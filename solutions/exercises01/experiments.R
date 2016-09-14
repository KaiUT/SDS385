"""

"""

source('~/Box Sync/PhDCourses/SDS385Statistical_models_for_big_data/SDS385/solutions/exercises01/Newton_method.R')

setwd('~/Box Sync/PhDCourses/SDS385Statistical_models_for_big_data/SDS385/data')
# import data
data <- read.csv('wdbc.csv', header=F)
y.BM <- data[ ,2]
y <- rep(0, length(y.BM))
y[y.BM == "B"] <- 1
y <- as.matrix(y)
x.variables <- as.matrix(data[ , 3:12])
x <- cbind(x.variables, rep(1,dim(x.variables)[1]))
colnames(x) <- NULL

# Gradient descent



# Newton's method
beta = matrix(0.0001, 11, 1)
betas <- my.newton(x, y, f.prime, f.prime2, beta, max.iter=50)
# plot

