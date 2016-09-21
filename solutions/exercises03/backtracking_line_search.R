source("~/Box Sync/PhDCourses/SDS385Statistical_models_for_big_data/SDS385/solutions/exercises01/gradient_descent.R")  # call functions `f.prime`, `l.beta` from `gradient_descent.R`.


# function for B
BFGS = function(x, y, betas0, betas.new, H) {
    delta = f.prime(x, y, betas.new) - f.prime(x, y, betas)
    s = betas.new - betas  # alpha * p.k
    gamma = 1 / (t(delta) %*% s)
    H.new = (diag(dim(betas0)[1]) - gamma[1,1] * s %*% t(delta)) %*% H %*%
        (diag(dim(betas0)[1]) - gamma[1,1] * delta %*% t(s)) + gamma[1,1] * s %*% t(s)
    return (H.new)
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
    l.prime = f.prime(x, y, betas)
    p.k = - H %*% l.prime
    return (p.k)
}


backtracking = function(x, y, betas0, alpha0, rho0, c0, H, BFGS=F) {
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
    l.prime = f.prime(x, y, betas0)
    l = l.beta(x, y, betas0)  # likelihood
    betas.new = betas0 + alpha * p.k
    l.new = l.beta(x, y, betas.new)
    if (BFGS) {
        H.new = BFGS(x, y, betas0, betas.new, H)
        H = H.new
        p.k = searchDirection(x, y, betas0, H)  # search direction
    }
    while (l.new > l + c * alpha * t(l.prime) %*% p.k) {
        alpha = rho * alpha
        betas.new = betas0 + alpha * p.k
        l.new = l.beta(x, y, betas.new)
    }
    step.size = alpha
    if (BFGS) {
        return (list(step.size, p.k))
    }
    return (step.size)
}


GDBacktracking = function(x, y, beta0, max.iter=50, alpha0, rho0, c0, H) {
    iterations = 0
    beta = beta0
    beta.matrix = beta
    l.tracking = c()
    while (iterations < max.iter) {
        iterations = iterations + 1
        step.size = backtracking(x, y, beta, alpha0, rho0, c0, H, BFGS=F)
        new.beta = beta - step.size * f.prime(x, y, beta)
        beta.matrix = cbind(beta.matrix, new.beta)
        l = l.beta(x, y, new.beta)
        l.tracking = append(l.tracking, l)
        beta = new.beta
    }
    return (list(beta.matrix, l.tracking))
}


NewtonBacktracking = function(x, y, beta0, max.iter=50, alpha0, rho0, c0, H) {
    iterations = 0
    beta = beta0
    beta.matrix = beta
    l.tracking = c()
    while (iterations < max.iter) {
        iterations = iterations + 1
        results = backtracking(x, y, beta, alpha0, rho0, c0, H, BFGS=T)
        step.size = results[[1]]
        p.k = results[[2]]
        new.beta = beta + step.size * p.k
        beta.matrix = cbind(beta.matrix, new.beta)
        l = l.beta(x, y, new.beta)
        l.tracking = append(l.tracking, l)
        beta = new.beta
    }
    return (list(beta.matrix, l.tracking))
}
