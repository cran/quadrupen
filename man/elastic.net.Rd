\name{elastic.net}
\alias{elastic.net}
\title{Fit a linear model with elastic-net regularization}
\usage{
  elastic.net(x, y, lambda1 = NULL, lambda2 = 0.01,
    penscale = rep(1, p), weights = rep(1, n),
    struct = 1 * bandSparse(p, k = 0), intercept = TRUE,
    normalize = TRUE, naive = FALSE,
    nlambda1 = ifelse(is.null(lambda1), 100, length(lambda1)),
    min.ratio = ifelse(n < p, 0.01, 0.005),
    max.feat = as.integer(ifelse(lambda2 < 0.01, min(n, p), min(3 * n, p))),
    control = list(), checkargs = TRUE)
}
\arguments{
  \item{x}{matrix of features, possibly sparsely encoded
  (experimental). Do NOT include intercept.}

  \item{y}{response vector.}

  \item{lambda1}{sequence of decreasing l1-penalty levels.
  If \code{NULL} (the default), a vector is generated with
  \code{nlambda1} entries, starting from a guessed level
  \code{lambda1.max} where only the intercept is included,
  then shrunken to \code{min.ratio*lambda1.max}.}

  \item{lambda2}{real scalar; tunes the ridges penalty in
  the Elastic-net. Default is 0.01. Set to 0 to recover the
  Lasso.}

  \item{penscale}{vector with real positive values that
  weight the l1-penalty of each features. Default set all
  weights to 1.}

  \item{struct}{matrix structuring the coefficients,
  possibly sparsely encoded (MUST be positive definite).
  Default is the identity matrix. See details below.}

  \item{weights}{observation weights for weighted
  least-squares. Default is 1 for each observation.}

  \item{intercept}{logical; indicates if an intercept
  should be included in the model. Default is \code{TRUE}.}

  \item{normalize}{logical; indicates if variables should
  be normalized to have unit L2 norm before fitting.
  Coefficients will then be rescaled to the original scale.
  Default is \code{TRUE}.}

  \item{naive}{logical; Compute either 'naive' of classic
  elastic-net as defined in Zou and Hastie (2006). Default
  is FALSE. See details below.}

  \item{nlambda1}{integer that indicates the number of
  values to put in the \code{lambda1} vector.  Ignored if
  \code{lambda1} is provided.}

  \item{min.ratio}{minimal value of l1-penalty that will be
  tried, as a fraction of the maximal \code{lambda1} value.
  A too small value of \code{lambda1} might lead to
  unstable path of solution when n<p.  The default value
  tries to avoid this, adapting to the 'n<p' context.
  Ignored if \code{lambda1} is provided.}

  \item{max.feat}{integer; limits the number of features
  ever to enter the model: the algorithm stops if this
  number is exceeded and \code{lambda1} is cut at the
  corresponding level. Default is
  \code{min(nrow(x),ncol(x))} for small \code{lambda2}
  (<0.01) and \code{min(3*nrow(x),ncol(x))} otherwise. Use
  with care, as it considerably changes the computation
  time.}

  \item{control}{list of argument controlling low level
  options of the algorithm --use with care and at your own
  risk-- : \itemize{%

  \item{\code{verbose}:}{logical; activate verbose mode.
  Default is \code{FALSE}.} \item{\code{timer}:}{logical;
  use to record the timing of the algorithm. Default is
  \code{FALSE}.} \item{\code{zero}:}{a practical zero.
  Default is \code{.Machine$double.eps}}
  \item{\code{max.iter}:}{the maximal number of iteration
  used to solve the problem for a given value of lambda1.
  Default is 500.} \item{\code{method}:}{a string for the
  underlying solver used. Either \code{"quadra"},
  \code{"pathwise"} or \code{"fista"}. Default is
  \code{"quadra"}} \item{\code{threshold}:}{a threshold for
  convergence. The algorithm stops when the optimality
  conditions are fulfill up to this threshold. Default is
  \code{1e-7} for \code{"quadra"} and \code{1e-2} for the
  first order methods.} \item{\code{monitor}:}{indicates if
  a monitoring of the convergence should be recorded, by
  computing a lower bound between the current solution and
  the optimum: when \code{'0'} (the default), no monitoring
  is provided; when \code{'1'}, the bound derived in
  Grandvalet et al. is computed; when \code{'>1'}, the
  Fenchel duality gap is computed along the algorithm.}}}

  \item{checkargs}{logical; should arguments be checked to
  (hopefully) avoid internal crashes? Default is
  \code{TRUE}.}
}
\value{
  a object with class \code{quadrupen}, see the
  documentation page \code{\linkS4class{quadrupen}} for
  details.
}
\description{
  Adjust a linear model with elastic-net regularization,
  mixing a (possibly weighted) l1-penalty (LASSO) and a
  (possibly structured) l2-penalty (ridge-like). The
  solution path is computed at a grid of values for the
  l1-penalty, fixing the value of the l2-penalty. See
  details for the criteria optimized.
}
\note{
  The optimized criterion is the following:
  \if{latex}{\deqn{% \hat{\beta}_{\lambda_1,\lambda_1} =
  \arg \min_{\beta} \frac{1}{2} (y - X \beta)^T D (y - X
  \beta) + \lambda_1 \|w \circ \beta \|_1 +
  \frac{\lambda_2}{2} \beta^T S \beta, }} \if{html}{\out{
  <center> &beta;<sup>hat</sup>
  <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
  argmin<sub>&beta;</sub> 1/2 RSS(&beta;,weights) +
  &lambda;<sub>1</sub> &#124; &beta; &#124;<sub>1</sub> +
  &lambda;/2 <sub>2</sub> &beta;<sup>T</sup> S &beta;,
  </center> }} \if{text}{\deqn{beta.hat(lambda1, lambda2) =
  argmin_beta 1/2 RSS(beta,weights) + lambda1 ||w.beta||1 +
  lambda2 beta' S beta,}} where the vector w is the
  \code{penscale} argument and S is the \code{struct}
  argument, a positive definite matrix (possibly of class
  \code{Matrix}). \if{latex}{\eqn{D} is a diagonal matrix
  whose diagonal corresponds to the 'weights' vector to
  perform weighted least-squares.}

  If \code{naive} is \code{TRUE}, the coefficients are not
  rescale. When \code{naive} is \code{FALSE} (the default)
  their are rescale by a scalar \code{1+lambda2}, as
  suggested in Zou and Hastie, to obtained the so-called
  Elastic-net estimator.
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

## This gives a great advantage to the elastic-net
## for support recovery
beta.lasso <- slot(crossval(x,y,lambda2=0) , "beta.min")
beta.enet  <- slot(crossval(x,y,lambda2=10), "beta.min")

cat("\\nFalse positives for the Lasso:", sum(sign(beta) != sign(beta.lasso)))
cat("\\nFalse positives for the Elastic-net:", sum(sign(beta) != sign(beta.enet)))
cat("\\nDONE.\\n")
}
\seealso{
  See also \code{\linkS4class{quadrupen}},
  \code{\link{plot.quadrupen}} and \code{\link{crossval}}.
}
\keyword{models,}
\keyword{regression}

