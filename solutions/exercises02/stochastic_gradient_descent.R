w.prob = function(x, beta) {
    sample.size = dim(x)[1]
    w = c()
    for (i in 1:sample.size){
        value = 1 / (1 + exp(-x[i, ] %*% beta))
        w = append(w, value)
    }
    return (w)
}


# tracking l(beta)
l.beta <- function(x, y, beta) {
    w = w.prob(x, beta)
    l = 0
    for (i in 1:length(y)) {
        value = y[i, ] * log((w[i] + 0.001) / (1 - w[i] + 0.001)) + 1 * log(1 - w[i] + 0.001)
        l = l + value
    }
    return (-l)
}


# tracking moving average l(beta)_t
l.beta.single <- function(x.sample, y.sample, beta) {
    w = w.prob(x.sample, beta)
    value = y.sample * log((w + 0.001) / (1 - w + 0.001)) + 1 * log(1 - w + 0.001)
    return (-value)
}


# sample data
sample.data <- function(x, y, max.iter, replace=T) {
    if (!replace & max.iter > length(y)) {
        a = max.iter %/% length(y)
        b = max.iter %% length(y)
        index = c()
        for (i in 1:a) {
            index1 = sample(1:length(y), length(y), replace=replace)
            index = c(index, index1)
        }
        index1 = sample(1:length(y), b, replace=replace)
        index = c(index, index1)
    }
    else {
        index = sample(1:length(y), max.iter, replace=replace)
    }
    return (index)
}


# SGD with constant step size.
SGD.constant.stepsize <- function(x, y, beta0, step.size, max.iter, replace=T, lambda=0.4) {
    count = 0
    l.total.tracking = c()
    l.average.tracking = c()
    l.weighted.tracking = c()
    l.average = 0
    l.sum = 0
    beta = beta0
    beta.matrix = beta
    # sample data
    index = sample.data(x, y, max.iter, replace=T)
    # evaluate betas
    for (i in index) {
        count = count + 1
        x.sample = x[i, , drop=F]
        y.sample = y[i, , drop=F]
        w = w.prob(x.sample, beta)
        #step.size is constant here.
        new.beta = beta - step.size * t(x.sample) %*% (1 * w - y.sample)
        # tracking total l(beta).
        l.total = l.beta(x, y, new.beta)
        l.total.tracking = append(l.total.tracking, l.total)
        # tracking exponentially weighted average
        l.single = l.beta.single(x.sample, y.sample, new.beta)
        l.weighed = lambda * l.single + (1 - lambda) * l.sum
        l.weighted.tracking = append(l.weighted.tracking, l.weighed)
        if (count > 1000 & count <= max.iter) {
            l.sum = l.sum + l.single
        }
        # tracking average l(beta)
        # l.average = l.sum / count
        l.average = (l.single + 2*l.average)/count
        l.average.tracking = append(l.average.tracking, l.average)

        beta.matrix = cbind(beta.matrix, new.beta)
        beta = new.beta
    }
    return (list(beta.matrix, l.total.tracking, l.weighted.tracking, l.average.tracking))
}


# Robbins-Monro rule.
RM.stepsize <- function(t, C, t0=1, alpha) {
    step.size = C * (t + t0) ** (-alpha)
    return (step.size)
}


# SGD with decaying step size.
SGD.RM.stepsize <- function(x, y, beta0, C, t0, alpha, max.iter, replace=T, lambda=0.8) {
    count = 0
    l.total.tracking = c()
    l.average.tracking = c()
    l.weighted.tracking = c()
    l.sum = 0
    beta = beta0
    beta.matrix = beta
    # sample data
    index = sample.data(x, y, max.iter, replace=T)
    # evaluate betas
    for (i in index) {
        count = count + 1
        x.sample = x[i, , drop=F]
        y.sample = y[i, , drop=F]
        w = w.prob(x.sample, beta)
        # RM rule for step size.
        step.size = RM.stepsize(count, C, t0, alpha)
        new.beta = beta - step.size * t(x.sample) %*% (1 * w - y.sample)
        # tracking total l(beta).
        l.total = l.beta(x, y, new.beta)
        l.total.tracking = append(l.total.tracking, l.total)
        # tracking exponentially weighted average
        l.single = l.beta.single(x.sample, y.sample, new.beta)
        l.weighed = lambda * l.single + (1 - lambda) * l.sum
        l.weighted.tracking = append(l.weighted.tracking, l.weighed)
        # tracking average l(beta)
        l.sum = l.sum + l.single
        l.average = l.sum / count
        l.average.tracking = append(l.average.tracking, l.average)

        beta.matrix = cbind(beta.matrix, new.beta)
        beta = new.beta
    }
    return (list(beta.matrix, l.total.tracking, l.weighted.tracking, l.average.tracking))
}


