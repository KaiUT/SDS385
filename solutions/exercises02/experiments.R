source('~/Box Sync/PhDCourses/SDS385Statistical_models_for_big_data/SDS385/solutions/exercises02/stochastic_gradient_descent.R')
source('~/Box Sync/PhDCourses/SDS385Statistical_models_for_big_data/SDS385/solutions/exercises01/Newton_method.R')

# import data
setwd('~/Box Sync/PhDCourses/SDS385Statistical_models_for_big_data/SDS385/data')
data <- read.csv('wdbc.csv', header=F)
y.BM <- data[ ,2]
y <- rep(0, length(y.BM))
y[y.BM == "B"] <- 1
y <- as.matrix(y)
x.variables <- as.matrix(data[ , 3:12])
x <- cbind(x.variables, rep(1,dim(x.variables)[1]))
colnames(x) <- NULL

# experiment for SGD with constant step size
betas <- matrix(0.0001, 11, 1)
betas.sgd.c <- SGD.constant.stepsize(x, y, betas, step.size=0.02, max.iter=1000, replace=T, lambda=0.8)

# experiment for SGD with constant step size
betas <- matrix(0.0001, 11, 1)
betas.sgd.RM <- SGD.RM.stepsize(x, y, betas, C=6, t0=2, alpha=1, max.iter=500, replace=T, lambda=0.8)

# Newton's method
betas <- matrix(0.0001, 11, 1)
betas.newton <- my.newton(x, y, f.prime, f.prime2, betas, max.iter=50)

# glm function in R
betas.glm <- glm(y[,1]~x[,1]+x[,2]+x[,3]+x[,4]+x[,5]+x[,6]+x[,7]+x[,8]+x[,9]+x[,10], family=binomial)



