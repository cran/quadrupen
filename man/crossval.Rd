\name{crossval}
\alias{crossval}
\title{Cross-validation function for quadrupen fitting methods.}
\usage{
  crossval(x, y, penalty = c("elastic.net", "bounded.reg"),
    K = 10,
    folds = split(sample(1:nrow(x)), rep(1:K, length = nrow(x))),
    lambda2 = 0.01, verbose = TRUE,
    mc.cores = detectCores(), ...)
}
\arguments{
  \item{penalty}{a string for the fitting procedure used
  for cross-validation. Either \code{"elastic.net"} or
  \code{"bounded.reg"}, at the moment. Default is
  \code{elastic.net}.}

  \item{x}{matrix of features, possibly sparsely encoded
  (experimental). Do NOT include intercept.}

  \item{y}{response vector.}

  \item{K}{integer indicating the number of folds. Default
  is 10.}

  \item{folds}{list of \code{K} vectors that describes the
  folds to use for the cross-validation. By default, the
  folds are randomly sampled with the specified K. The same
  folds are used for each values of \code{lambda2}.}

  \item{lambda2}{tunes the \eqn{\ell_2}{l2}-penalty
  (ridge-like) of the fit. If none is provided, the default
  scalar value of the corresponding fitting method is used
  and a simple CV is performed. If a vector of values is
  given, double cross-validation is performed (both on
  \code{lambda1} and \code{lambda2}, using the same folds
  for each \code{lambda2}).}

  \item{verbose}{logical; indicates if the progression (the
  current lambda2) should be displayed. Default is
  \code{TRUE}.}

  \item{mc.cores}{the number of cores to use. The default
  uses all the cores available.}

  \item{...}{additional parameters to overwrite the
  defaults of the fitting procedure identified by the
  \code{'penalty'} argument. See the corresponding
  documentation (\code{\link{elastic.net}} or
  \code{\link{bounded.reg}}).}
}
\value{
  An object of class "cvpen" for which a \code{plot} method
  is available.
}
\description{
  Function that computes K-fold (double) cross-validated
  error of a \code{quadrupen} fit. If no \code{lambda2} is
  provided, simple cross validation on the \code{lambda1}
  parameter is performed. If a vector \code{lambda2} is
  passed as an argument, double cross-validation is
  performed.
}
\note{
  If the user runs the fitting method with option
  \code{'bulletproof'} set to \code{FALSE}, the algorithm
  may stop at an early stage of the path. Early stops are
  handled internally, in order to provide results on the
  same grid of penalty tuned by \eqn{\lambda_1}{lambda1}.
  This is done by means of \code{NA} values, so as mean and
  standard error are consistently evaluated. If, while
  cross-validating, the procedure experiences too much
  early stoppings, a warning is sent to the user, in which
  case you should reconsider the grid of \code{lambda1}
  used for the cross-validation.  If \code{bulletproof} is
  \code{TRUE} (the default), there is nothing to worry
  about, except a possible slow down when any switching to
  the proximal algorithm is required.
}
\examples{
\dontrun{
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor  <- 0.75
Soo  <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variable
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
\seealso{
  \code{\linkS4class{quadrupen}}, \code{\link{plot.cvpen}}
  and \code{\linkS4class{cvpen}}.
}
\keyword{models,}
\keyword{regression}

