\name{crossval}
\alias{crossval}
\title{Cross-validation function for quadrupen fitting methods.}
\usage{
  crossval(x, y, method = "elastic.net", K = 10,
    folds = split(sample(1:nrow(x)), rep(1:K, length = nrow(x))),
    lambda2 = 0.01, verbose = TRUE, ...)
}
\arguments{
  \item{method}{a string for the fitting procedure used for
  cross-validation. Only \code{"elastic.net"} is available
  for the moment}

  \item{x}{matrix of features, possibly sparsely encoded
  (experimental). Do NOT include intercept.}

  \item{y}{response vector.}

  \item{K}{integer indicating the number of folds. Default
  is 10.}

  \item{folds}{list of \code{K} vectors that describes the
  folds to use for the cross-validation. By default, the
  folds are randomly sampled with the specified K. The same
  folds are used for each values of \code{lamdba2}.}

  \item{lambda2}{tunes the l2-penalty (ridge-like) of the
  fit. If none is provided, the default scalar value of the
  corresponding fitting method is used and a simple CV is
  performed. If a vector of values is provided, double
  cross-validation is performed (both on lambda1 and
  lambda2, using the same folds for each lambda2).}

  \item{verbose}{logical; indicates if the progression (the
  current fold number) should be displayed. Default is
  \code{TRUE}.}

  \item{...}{additional parameters to overwrite the
  defaults of the \code{'method'} fitting procedure. See
  the corresponding documentation (e.g.,
  \code{\link{elastic.net}})}
}
\value{
  An object of class "cvpen" for which a \code{plot} method
  is available.
}
\description{
  Function that computes K-fold (double) cross-validated
  error of a \code{quadrupen} fit. If no lambda2 is
  provided, simple cross validation on the lambda1
  parameter is performed. If a vector lambda2 is provided,
  double cross-validation is performed.
}
\examples{
rm(list=ls())
library(quadrupen)
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) \%*\% chol(Sigma))
y <- 10 + x \%*\% beta + rnorm(n,0,10)

## Use fewer lambda1 values by overwritting the default parameters
## and cross-validate over the sequences lambda1 and lambda2
cvout <- crossval(x,y, lambda2=10^seq(2,-2,len=50), nlambda1=50)

beta.1se <- slot(cvout, "beta.1se")
beta.min <- slot(cvout, "beta.min")

cat("\\nFalse positives with the minimal CV choice: ", sum(sign(beta) != sign(beta.min)))
cat("\\nFalse positives with the 1 standard-error rule: ", sum(sign(beta) != sign(beta.1se)))
}
\seealso{
  \code{\linkS4class{quadrupen}}, \code{\link{plot.cvpen}}
  and \code{\linkS4class{cvpen}}.
}
\keyword{models,}
\keyword{regression}

