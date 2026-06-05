rm(list=ls())
library(quadrupen)
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,5)

## Build a vector of label for true nonzeros
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- c("relevant")
labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))

## Call to stability selection function, 200 subsampling
enet <- elastic_net(x, y, lambda2 = 1e-2, minratio = 1e-3)
stab <- stability(enet, n_subsamples = 200)

## Build the plot an recover the selected variable for a given cutoff
## and per-family error rate
plot(stab, labels=labels)
stab$plot(labels=labels)
stab$plot(labels=labels, nvarsel=5)
stab$plot(sel_mode = "PFER", labels=labels)
stab$plot(xvar="fraction", sel_mode = "PFER", labels=labels, PFER = 2)
plot(stab, xvar="fraction", sel_mode = "PFER", labels=labels, PFER = 4)

cat("\n\nFalse positives for the randomized Elastic-net with stability selection: ",
     sum(labels[stab$selection()] != "relevant"))
cat("\nDONE.\n")

