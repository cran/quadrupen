library(quadrupen)
require(Matrix)

beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)
p <- ncol(x)
plot(beta, type = "s")

## Structuring matrix :
# Laplacian of a community graph
# (clique prior)
Ioo <- matrix(1,25,25)
Iww <- matrix(1,10,10)
A <- bdiag(Ioo,Iww,Ioo,Iww,Ioo)
diag(A) <- 0
L2 <- -A
diag(L2) <- colSums(A) + 1e-2

labels <- paste("segment", rep(1:5, c(25,10,25,10,25)), sep="-" )
lasso_fit <- lasso(x,y,intercept=F)
plot(lasso_fit, labels=labels)
fused_lasso_fit <- fused_lasso(x,y,lambda2=5,intercept=F)
plot(fused_lasso_fit, labels=labels)
ridge_fit <- ridge(x,y,intercept=F)
plot(ridge_fit, labels=labels)
ridge_struct <- ridge(x,y, struct=L2, lambda_max=1000)
plot(ridge_struct, labels=labels)
enet_fit <- elastic_net(x,y, lambda2 = 5, intercept=F)
plot(enet_fit, labels=labels)
enet_struct <- elastic_net(x,y, lambda2 = 5, struct=L2)
plot(enet_struct, labels=labels)
enet_struct$debias <- TRUE
plot(enet_struct, labels=labels)

beta_lasso <- as.numeric(lasso(x,y,lambda=10)$coefficients)
beta_fused_lasso <- as.numeric(fused_lasso(x,y,lambda1=20, lambda2 = 5)$coefficients)
beta_ridge <- as.numeric(ridge(x,y,lambda=1)$coefficients)
beta_ridge_struct <- as.numeric(ridge(x,y,lambda=10,struct=L2)$coefficients)
beta_enet_struct  <- as.numeric(elastic_net(x,y,lambda2=5,lambda1=20,struct=L2)$coefficients)
beta_enet_struct_refit  <- as.numeric(elastic_net(x,y,lambda2=5,lambda1=20,struct=L2,refit=TRUE)$coefficients)

par(mfrow=c(3,2))
plot(beta_lasso, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_fused_lasso, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_ridge, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_ridge_struct, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_enet_struct, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_enet_struct_refit, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
