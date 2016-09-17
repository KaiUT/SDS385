library('Matrix')

# inversion method
inversion.method = function(X, y, W) {
    beta = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
    return (beta)
}


# Cholesky method
Cholesky.method = function(X, y, W) {
    A = t(X)%*%W%*%X
    b = t(X)%*%W%*%y
    L = Matrix::chol(A)
    beta = solve(L) %*% solve(t(L)) %*% b
    return(beta)
}
