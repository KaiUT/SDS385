lasso_loglik <- function (X, y, betas, lambda) {
    #
    #
    #
    l = 1/(2*nrow(X)) * t(y - X %*% betas) %*% (y - X %*% betas)
    phi = lambda * sum(abs(betas))
    return (l + phi)
}


lasso_penalty <- function(mu, gamma, lambda) {
    #
    #
    #
    betas.new = matrix(0, nrow=nrow(mu))
    for (feature in 1:nrow(mu)) {
        betas.new[feature, 1] = sign(mu[feature, 1]) * max(0, abs(mu[feature, 1]) - gamma * lambda)
    }
    return (betas.new)
}


proximal_gd <- function (X, y, betas0, lambda, gamma, maxiter) {
    # X, y should be matrix
    # X should be normalized.

    # initialize parameters
    betas = betas0
    lambda = lambda
    gamma = gamma
    nfeatures = ncol(X)
    nsamples = nrow(X)
    # track betas, and log likelihood values
    betas.matrix <- matrix(0, nrow=nfeatures, ncol=maxiter+1)
    loglik <- rep(0, maxiter+1)
    betas.matrix[, 1] <- betas
    loglik[1] <- lasso_loglik(X, y, betas, lambda)

    # for loop to update betas
    for (iter in 1:maxiter) {
        # calculate gradient descent
        l.prime = 1/(2*nsamples) * (-2 * t(y) %*% X + 2 * t(betas) %*% t(X) %*% X)
        # update mu
        mu = betas - gamma * t(l.prime)
        # update betas
        betas = lasso_penalty(mu, gamma, lambda)
        # track betas and log likelihood
        betas.matrix[, iter+1] = betas
        loglik[iter+1] = lasso_loglik(X, y, betas, lambda)
    }
    return (list(betas.matrix, loglik))
}



acc_proximal_gd <- function (X, y, betas0, lambda, gamma, maxiter) {
    # X, y should be matrix
    # X should be normalized.

    # initialize parameters
    betas = betas0
    lambda = lambda
    gamma = gamma
    nfeatures = ncol(X)
    nsamples = nrow(X)
    # track betas, and log likelihood values
    betas.matrix <- matrix(0, nrow=nfeatures, ncol=maxiter+1)
    loglik <- rep(0, maxiter+1)
    betas.matrix[, 1] <- betas
    loglik[1] <- lasso_loglik(X, y, betas, lambda)
    # store s and z
    s = rep(0, maxiter+1)
    z = matrix(0, nrow=nfeatures, ncol=maxiter+1)
    z[, 1] <- betas0

    # for loop to update betas. Start from iteration 2.
    for (iter in 1:maxiter) {
        # calculate gradient descent
        l.prime = 1/(2*nsamples) * (-2 * t(y) %*% X + 2 * t(betas.matrix[, iter]) %*% t(X) %*% X)
        # update mu
        mu = betas.matrix[, iter] - gamma * t(l.prime)
        # update betas
        z[, iter+1] = lasso_penalty(mu, gamma, lambda)
        s[iter+1] = (1 + sqrt(1 + 4 * s[iter]^2)) / 2
        betas = z[, iter+1] + (s[iter] - 1)/s[iter+1] *(z[, iter+1] - z[, iter])
        # track betas and log likelihood
        betas.matrix[, iter+1] = betas
        loglik[iter+1] = lasso_loglik(X, y, betas, lambda)
    }
    return (list(betas.matrix, loglik))
}
