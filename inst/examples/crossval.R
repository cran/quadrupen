rm(list=ls())
library(quadrupen)
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor  <- 0.75
Soo  <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variable
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.1
diag(Sigma) <- 1
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

  
## Use fewer lambda1 values by overwritting the default parameters
## and cross-validate over the sequences lambda1 and lambda2
enet <- elastic_net(x, y, nlambda1 = 40, minratio = 0.01)
cvGRID <- cross_validate(enet, lambda2=10^seq(.5,-1,len=40))
plot(cvGRID)

## Rerun simple cross-validation with the appropriate lambda2
enet <- elastic_net(x, y, lambda2 = cvGRID$lambda2_min)
cv10K <- cross_validate(enet)
plot(cv10K)
beta10K <- enet$get_model("CV_min")[-1]

## Try leave one out also
cvLOO <- cross_validate(enet, K=n)
plot(cvLOO)
betaLOO <- enet$get_model("CV_min")[-1]

## Performance for selection purpose
cat("\nFalse positives with the minimal 10-CV choice: ", sum(sign(beta) != sign(beta10K)))
cat("\nFalse positives with the minimal LOO-CV choice: ", sum(sign(beta) != sign(betaLOO)))
