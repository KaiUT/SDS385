w.prob = function(x, beta) {
    w = 1 / (1 + exp(-x %*% beta))
    return (w)
}


# tracking l(beta)
l.beta <- function(x, y, beta) {
    l = t(y) %*% x %*% beta - t(matrix(1, nrow=dim(y)[1], 1)) %*% log(1 + exp(x %*% beta))
    return (-l)
}


f.prime = function(x, y, beta){
    w = w.prob(x, beta)
    # w.matrix = as.matrix(w)
    value = -t(x) %*% (y - 1 * w)
    return (value)
}


gradient.descent = function(x, y, beta0, step.size= 0.1, max.iter=50) {
    iterations = 0
    beta = beta0
    beta.matrix = beta
    l.tracking = c()
    while (iterations < max.iter) {
        iterations = iterations + 1
        new.beta = beta - step.size * f.prime(x, y, beta)
        beta.matrix = cbind(beta.matrix, new.beta)
        l = l.beta(x, y, new.beta)
        l.tracking = append(l.tracking, l)
        beta = new.beta
    }
    return (list(beta.matrix, l.tracking))
}
