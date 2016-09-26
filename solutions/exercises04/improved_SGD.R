# Improved stochastic gradient descent method.

library("Matrix")


wProb <- function(x, betas) {
    w = 1 / (1 + exp(-x %*% betas))
    return (w)
}


# # tracking log likelihood
# logLik <- function(x, y, betas) {
    # l = t(y) %*% x %*% betas - t(matrix(1, nrow=dim(y)[1], 1)) %*% log(1 + exp(x %*% betas))
    # return (-l)
# }


# logLikPrime = function(x, y, betas){
    # w = wProb(x, betas)
    # # w.matrix = as.matrix(w)
    # value = -t(x) %*% (y - 1 * w)
    # return (value)
# }


# objective function with L2 penalty
logLik <- function(x, y, betas, lambda=0) {
    l = t(y) %*% x %*% betas - t(matrix(1, nrow=dim(y)[1], 1)) %*% log(1 + exp(x %*% betas)) + lambda * crossprod(betas)
    return (-l)
}


# log likelihood prime with L2 penalty
logLikPrime = function(x, y, betas, lambda=0){
    w = wProb(x, betas)
    # w.matrix = as.matrix(w)
    value = -t(x) %*% (y - 1 * w) + lambda * 2 * t(betas)
    return (value)
}


# sample data
samplingData <- function(x, y, max.iter, replace=T) {
    if (!replace & max.iter > length(y)) {
        a = max.iter %/% length(y)
        b = max.iter %% length(y)
        index = c()
        for (i in 1:a) {
            index.temp = sample(1:length(y), length(y), replace=replace)
            index = c(index, index.temp)
        }
        index.temp = sample(1:length(y), b, replace=replace)
        index = c(index, index.temp)
    } else {
        index = sample(1:length(y), max.iter, replace=replace)
    }
    return (index)
}


# sample a small subset of data as a calibration set.
minibatch <- function(x, y, batch.size) {
    index = sample(1:length(y), batch.size)
    y.subset = y[index, , drop=F]
    x.subset = x[index, , drop=F]
    return (list(x.subset, y.subset))
}


searchDirection = function(x, y, betas, H) {
    # Determine the search direction.
    #
    # Args:
    #   x:
    #   y:
    #   betas:
    #   B:
    #
    # Returns:
    #
    l.prime = logLikPrime(x, y, betas) / length(y)
    p.k = - H %*% l.prime
    return (p.k)
}


backtracking = function(x, y, betas0, H, alpha0=1, rho0=0.2, c0=1e-04, BFGS=F) {
    # Backtracking line search method
    #
    # Args:
    #   x:
    #   y:
    #   betas:
    #   B:
    #
    # Returns:
    #
    alpha = alpha0
    rho = rho0
    c= c0
    p.k = searchDirection(x, y, betas0, H)  # search direction
    l.prime = logLikPrime(x, y, betas0)
    l = logLik(x, y, betas0)  # likelihood
    betas.new = betas0 + alpha * p.k
    l.new = logLik(x, y, betas.new)
    # if (BFGS) {
        # H.new = BFGS(x, y, betas0, betas.new, H)
        # H = H.new
        # p.k = searchDirection(x, y, betas0, H)  # search direction
    # }
    while (l.new > l + c * alpha * t(l.prime) %*% p.k) {
        alpha = rho * alpha
        betas.new = betas0 + alpha * p.k
        l.new = logLik(x, y, betas.new)
    }
    step.size = alpha
    # if (BFGS) {
        # return (list(step.size, p.k))
    # }
    return (step.size)
}


# exponentially moving weighted average
EMWA <- function(value.current, avg.prev, weight=0.5) {
    avg = weight* value.current + (1 - weight) * avg.prev
    return (avg)
}


LineSearchSGD <- function(x, y, betas0, H, max.iter, replace=T, refresh.rate=10, batch.size=10, weight=0.5, ...) {
    # Stochastic gradient descent method with changing step size.
    # Args:

    # sample data
    index = samplingData(x, y, max.iter, replace=replace)
    x.sample = x[index, , drop=F]
    y.sample = y[index, , drop=F]
    # initialization
    betas = betas0
    betas.matrix = matrix(0, nrow=length(betas), ncol=length(index)+1)
    betas.matrix[, 1] = betas
    l.average = 0
    l.weighted = 0
    l.average.tracking = rep(0, length(index))
    l.weighted.tracking = rep(0, length(index))

    # evaluate betas
    for (i in 1:length(index)) {
        x.single = x.sample[i, , drop=F]
        y.single = y.sample[i, , drop=F]
        # sample a minibatch to determine next step size.
        remainder = (i %% refresh.rate == 1 | i / refresh.rate == i)
        if (remainder) {
            data.minibatch = minibatch(x.sample, y.sample, batch.size)
            step.size = backtracking(data.minibatch[[1]], data.minibatch[[2]], betas0, H, ...)
        }
        # calculate new betas
        new.betas = betas - step.size * logLikPrime(x.single, y.single, betas)
        # calculate log likelihood for a single data point.
        l.single = logLik(x.single, y.single, new.betas)
        # tracking exponentially weighted average
        if (i > 1000 & i <= max.iter) {
            l.weighted = EMWA(l.single, l.weighted, weight=weight)
        } else {
            l.weighted = l.single
        }
        l.weighted.tracking[i] = l.weighted
        # tracking average l(beta)
        l.average = (l.single + (i - 1) * l.average) / i
        l.average.tracking[i] = l.average
        # add updated betas to matrix
        betas.matrix[, i+1] = new.betas
        betas = new.betas
    }
    return (list(betas.matrix,l.weighted.tracking, l.average.tracking))
}


AdaGradSGD <- function(x, y, betas0, step.size, max.iter, replace=T, weight=0.5) {
    # Adaptive gradient stochastic gradient descent.
    # Args:

    # sample data
    index = samplingData(x, y, max.iter, replace=replace)
    x.sample = x[index, , drop=F]
    y.sample = y[index, , drop=F]
    # initialization
    betas = betas0
    betas.matrix = matrix(0, nrow=length(betas), ncol=length(index)+1)
    betas.matrix[, 1] = betas
    gd.cum = matrix(0, nrow=length(betas))
    l.average = 0
    l.weighted = 0
    l.average.tracking = rep(0, length(index))
    l.weighted.tracking = rep(0, length(index))

    # evaluate betas
    for (i in 1:length(index)) {
        x.single = x.sample[i, , drop=F]
        y.single = y.sample[i, , drop=F]
        # adjust step size
        gd.single = logLikPrime(x.single, y.single, betas)
        gd.cum = gd.cum + as.matrix(diag(tcrossprod(gd.single)))
        # gd.cum = gd.cum + gd.single * gd.single
        adjusted.gd = gd.single / (sqrt(gd.cum) + 1e-06)
        # calculate new betas
        new.betas = betas - step.size * adjusted.gd
        # calculate log likelihood for a single data point.
        l.single = logLik(x.single, y.single, new.betas)
        # tracking exponentially weighted average
        if (i > 1000 & i <= max.iter) {
            l.weighted = EMWA(l.single, l.weighted, weight=weight)
        } else {
            l.weighted = l.single
        }
        l.weighted.tracking[i] = l.weighted
        # tracking average l(beta)
        l.average = (l.single + (i - 1) * l.average) / i
        l.average.tracking[i] = l.average
        # add updated betas to matrix
        betas.matrix[, i+1] = new.betas
        betas = new.betas
    }
    return (list(betas.matrix,l.weighted.tracking, l.average.tracking))
}


# Create sparse matrix
sMatrix <- function(x) {
    x.length = length(x)
    i = sample(1:x.length, x.length)
    j = sample(1:x.length, x.length)
    sparse.matrix <- sparseMatrix(i=i, j=j, x=as.vector(x))
    return (sparse.matrix)
}
