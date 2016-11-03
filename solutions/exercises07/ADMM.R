objective_fn <- function(X, y, betas, lambda) {
    # X,y: scaled matrix
    # betas: vector

    value = 1/2 * sum((X %*% betas - y)^2) + lambda * sum(abs(betas))
    return (value)
}


soft_thresholding <- function(betas, u, lambda, rho) {
    # betas: vector
    # u: vector
    # z: vector

    value = betas + u
    z = sign(value) * pmax(betas + u - lambda/rho, rep(0, length(betas)))
    return (z)
}


ADMMLasso <- function(X, y, rho, lambda, maxiter) {
    # X, y must be scaled matrix
    # rho: step size
    # lambda: penalty
    # maxiter: number of iterations

    nsamples = nrow(X)
    nfeatures = ncol(X)
    lambda = nsamples * lambda  # rescale lambda to compare with function `glmnet`.

    # pre caching
    tXXrho = solve(crossprod(X) + diag(rep(rho, nfeatures)))
    Xty = crossprod(X, y)

    # initialize betas, zs and us
    betas.matrix = matrix(0, nrow=nfeatures, ncol=maxiter)
    z.matrix = matrix(0, nrow=nfeatures, ncol=maxiter+1)
    u.matrix = matrix(0, nrow=nfeatures, ncol=maxiter+1)
    objectives = rep(0, maxiter)

    for (i in 1:maxiter) {
        # update betas
        betas.matrix[, i] = tXXrho %*% (Xty + rho * (z.matrix[, i] - u.matrix[, i]))
        # update z
        z.matrix[, i+1] = soft_thresholding(betas.matrix[, i], u.matrix[, i], lambda, rho)
        # update u
        u.matrix[, i+1] = u.matrix[, i] + betas.matrix[, i] - z.matrix[, i+1]
        # update objective value
        objectives[i] = objective_fn(X, y, betas.matrix[, i], lambda)
    }
    return (list(betas=betas.matrix, betas.final=betas.matrix[, maxiter], objective = objectives))
}
