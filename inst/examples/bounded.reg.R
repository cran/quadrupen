rm(list=ls())
library(quadrupen)
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variable
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

## Infinity norm without/with an additional l2 regularization term
## and with structuiring prior
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- "relevant"

BREG0    <- bounded_reg(x,y,lambda2=0) 
BREGNET  <- bounded_reg(x,y,lambda2=10)
BREGSNET <- bounded_reg(x,y,lambda2=10 , struct=solve(Sigma))

plot(BREG0, label=labels) ## a mess
plot(BREGNET, label=labels) ## good guy are the boundaries
plot(BREGSNET, label=labels) ## even better

cv_breg0     <- cross_validate(BREG0)
cv_bregnet   <- cross_validate(BREGNET)
cv_bregnsnet <- cross_validate(BREGSNET)

plot(cv_breg0)
plot(cv_bregnet)
plot(cv_bregnsnet)

beta_breg0    <- BREG0$get_model("CV_1se")[-1]
beta_bregnet  <- BREGNET$get_model("CV_1se")[-1]
beta_bregsnet <- BREGSNET$get_model("CV_1se")[-1]




