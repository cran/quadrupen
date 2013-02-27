\docType{methods}
\name{plot,cvpen-method}
\alias{plot.cvpen}
\alias{plot,cvpen-method}
\title{Plot method for cross validated error of a \code{quadrupen} model}
\arguments{
  \item{x}{output of a \code{crossval} run (must be of
  class \code{cvpen}).}

  \item{y}{used for S4 compatibility.}

  \item{main}{the main title, with a hopefully appropriate
  default definition.}

  \item{xvar}{variable to plot on the X-axis: either
  \code{"lambda"} (\eqn{\lambda_1}{lambda1} penalty level)
  or \code{"fraction"} (fraction of \eqn{\ell_1}{l1}-norm
  of the shruken coefficients). Default is set to
  \code{"lambda"}.}

  \item{log.scale}{logical; indicates if a log-scale should
  be used when \code{xvar="lambda"}. Ignored for 2D
  cross-validation plot.}

  \item{plot}{logical; indicates if the graph should be
  plotted. Default is \code{TRUE}.}

  \item{reverse}{logical; should the X-axis by reversed
  when \code{xvar=lambda}? Default is FALSE.  Ignored for
  2D cross-validation plot.}

  \item{...}{used for S4 compatibility.}
}
\value{
  a \pkg{ggplot2} object which can be plotted via the
  \code{print} method.
}
\description{
  Produce a plot of the cross validated error of a
  \code{quadrupen} model.
}
\examples{
\dontrun{
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor  <- 0.75
Soo  <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.1
diag(Sigma) <- 1
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) \%*\% chol(Sigma))
y <- 10 + x \%*\% beta + rnorm(n,0,10)

## Use fewer lambda1 values by overwritting the default parameters
## and cross-validate over the sequences lambda1 and lambda2
cv.double <- crossval(x,y, lambda2=10^seq(2,-2,len=50), nlambda1=50)
## Rerun simple cross-validation with the appropriate lambda2
cv.10K <- crossval(x,y, lambda2=slot(cv.double, "lambda2.min"))
## Try leave one out also
cv.loo <- crossval(x,y, K=n, lambda2=slot(cv.double, "lambda2.min"))

plot(cv.double)
plot(cv.10K)
plot(cv.loo)

## Performance for selection purpose
beta.min.10K <- slot(cv.10K, "beta.min")
beta.min.loo <- slot(cv.loo, "beta.min")

cat("\\nFalse positives with the minimal 10-CV choice: ", sum(sign(beta) != sign(beta.min.10K)))
cat("\\nFalse positives with the minimal LOO-CV choice: ", sum(sign(beta) != sign(beta.min.loo)))
}
}

