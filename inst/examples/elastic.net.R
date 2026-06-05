rm(list=ls())
library(quadrupen)
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- "relevant"
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
n <- 100
p <- length(beta)
x <- as.matrix(matrix(rnorm(p*n),n,p) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

## Comparing the solution path of the LASSO, the Elastic-net and the
## Structured Elastic.net
LASSO <- lasso(x,y, minratio = 1e-2)
ENET  <- elastic_net(x,y, minratio = 1e-2, lambda2 = 2)
SENET <- elastic_net(x,y, minratio = 1e-2, lambda2 = 1, struct=solve(cov(x)) + diag(1e-3, p, p))

plot(LASSO, label=labels) ## a mess
plot(ENET , label=labels) ## a lot better
plot(SENET, label=labels) ## even better

## This gives a great advantage to the elastic-net
## for support recovery
cv_lasso <- cross_validate(LASSO)
cv_enet  <- cross_validate(ENET)
cv_senet <- cross_validate(SENET)

plot(cv_lasso)
plot(cv_enet)
plot(cv_senet)

beta_lasso <- LASSO$get_model("CV_1se")[-1]
beta_enet  <- ENET$get_model("CV_1se")[-1]
beta_senet <- SENET$get_model("CV_1se")[-1]

cat("\nFalse positives for the Lasso:", sum(sign(beta) != sign(beta_lasso)))
cat("\nFalse positives for the Elastic-net:", sum(sign(beta) != sign(beta_enet)))
cat("\nFalse positives for Structured Elastic-net:", sum(sign(beta) != sign(beta_senet)))
cat("\nDONE.\n")

## Now after debiasing
LASSO$debias <- TRUE
ENET$debias <- TRUE
SENET$debias <- TRUE
cv_lasso <- cross_validate(LASSO)
cv_enet  <- cross_validate(ENET)
cv_senet <- cross_validate(SENET)

plot(cv_lasso)
plot(cv_enet)
plot(cv_senet)

beta_lasso <- LASSO$get_model("CV_min")[-1]
beta_enet  <- ENET$get_model("CV_min")[-1]
beta_senet <- SENET$get_model("CV_min")[-1]

cat("\nFalse positives for the Lasso:", sum(sign(beta) != sign(beta_lasso)))
cat("\nFalse positives for the Elastic-net:", sum(sign(beta) != sign(beta_enet)))
cat("\nFalse positives for Structured Elastic-net:", sum(sign(beta) != sign(beta_senet)))
cat("\nDONE.\n")

