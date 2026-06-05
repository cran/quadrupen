## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

## Plot penalized criteria for the Elastic-net path
enet <- elastic_net(x,y, lambda2=1, refit = TRUE)
enet$information_criteria$plot()
enet$information_criteria$plot("BIC")
enet$information_criteria$plot(c("AIC", "BIC"))

## Plot penalized criteria for the Elastic-net path
enet <- elastic_net(x,y, lambda2=1, refit = FALSE)
plot(enet, "criteria")

## Plot penalized criteria for the Bounded regression
breg <- bounded_reg(x,y, lambda2=1)
plot(breg, "criteria")

## Plot penalized criteria for the Ridge regression
ridge <- ridge(x,y)
plot(ridge, "criteria")

## Plot penalized criteria for the LAVA path
lav <- lava(x, y, lambda2 = 0.1, minratio = )
lav$information_criteria$plot(c("AIC", "BIC", "mBIC"))
