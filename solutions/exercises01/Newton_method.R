w.prob = function(x, beta) {
    sample.size = dim(x)[1]
    w = c()
    for (i in 1:sample.size){
        value = 1 / (1 + exp(-x[i, ] %*% beta))
        w = append(w, value)
    }
    return (w)
}


f.prime = function(x, y, beta){
    w = w.prob(x, beta)
    w.matrix = as.matrix(w)
    value = -t(x) %*%(y - 1 * w)
    return (value)
}


f.prime2 = function(x, y, beta){
    w = w.prob(x, beta)
    W = diag(1 * w * (1 - w))
    value = as.matrix(t(x) %*% W %*% x)
    return (value)
}


my.newton = function(x, y, f.prime, f.prime2, beta0, max.iter=50) {
    iterations = 0
    beta = beta0
    beta.matrix = beta
    while (iterations < max.iter) {
        iterations = iterations + 1
        new.beta = beta - solve(f.prime2(x, y, beta)) %*% f.prime(x, y, beta)  # using default step size 1.
        beta.matrix = cbind(beta.matrix, new.beta)
        beta = new.beta
    }
    return (beta.matrix)
}
