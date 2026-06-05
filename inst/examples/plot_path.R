rm(list=ls())
library(quadrupen)
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
## This gives a great advantage to the elastic-net
## for support recovery; have a look at image(Sigma)
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)
lasso <- elastic_net(x,y, lambda2=0)
enet  <- elastic_net(x,y, lambda2=10)
breg1  <- bounded_reg(x,y, lambda2=1e-2)
breg2  <- bounded_reg(x,y, lambda2=10)

## Plot the Lasso path
lasso$plot_path(title="Lasso solution path")
## Plot the Elastic-net path
enet$plot_path(title = "Elastic-net solution path")
## Plot the Elastic-net path (fraction on X-axis, unstandardized coefficient)
enet$plot_path(standardize=FALSE, xvar="fraction")
## Plot the bounded regression path (fraction on X-axis)
breg1$plot_path()
## Plot another bounded regression path (fraction on X-axis, unstandardized coefficient)
breg2$plot_path(standardize=FALSE, xvar="fraction")
