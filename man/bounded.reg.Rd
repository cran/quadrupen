\name{bounded.reg}
\alias{bounded.reg}
\title{Fit a linear model with infinity-norm plus ridge-like regularization}
\usage{
  bounded.reg(x, y, lambda1 = NULL, lambda2 = 0.01,
    penscale = rep(1, p), struct = NULL, intercept = TRUE,
    normalize = TRUE, naive = FALSE,
    nlambda1 = ifelse(is.null(lambda1), 100, length(lambda1)),
    min.ratio = ifelse(n <= p, 0.01, 0.001),
    max.feat = ifelse(lambda2 < 0.01, min(n, p), min(4 * n, p)),
    control = list(), checkargs = TRUE)
}
\arguments{
  \item{x}{matrix of features, possibly sparsely encoded
  (experimental). Do NOT include intercept. When normalized
  os \code{TRUE}, coefficients will then be rescaled to the
  original scale.}

  \item{y}{response vector.}

  \item{lambda1}{sequence of decreasing
  \eqn{\ell_\infty}{l-infinity} penalty levels. If
  \code{NULL} (the default), a vector is generated with
  \code{nlambda1} entries, starting from a guessed level
  \code{lambda1.max} where only the intercept is included,
  then shrunken to \code{min.ratio*lambda1.max}.}

  \item{lambda2}{real scalar; tunes the
  \eqn{\ell_2}{l2}-penalty in the bounded regression.
  Default is 0.01. Set to 0 to regularize only by the
  infinity norm (be careful regarding numerical stability
  in that case, particularly in the high dimensional
  setting).}

  \item{penscale}{vector with real positive values that
  weight the infinity norm of each feature. Default set all
  weights to 1. See details below.}

  \item{struct}{matrix structuring the coefficients,
  possibly sparsely encoded.  Must be at least positive
  semidefinite (this is checked internally if the
  \code{checkarg} argument is \code{TRUE}).  If \code{NULL}
  (the default), the identity matrix is used. See details
  below.}

  \item{intercept}{logical; indicates if an intercept
  should be included in the model. Default is \code{TRUE}.}

  \item{normalize}{logical; indicates if variables should
  be normalized to have unit L2 norm before fitting.
  Default is \code{TRUE}.}

  \item{naive}{logical; Compute either 'naive' of 'classic'
  bounded regression: mimicing the Elastic-net, the vector
  of parameters is rescaled by a coefficient
  \code{(1+lambda2)} when \code{naive} equals \code{FALSE}.
  No rescaling otherwise. Default is \code{FALSE}.}

  \item{nlambda1}{integer that indicates the number of
  values to put in the \code{lambda1} vector.  Ignored if
  \code{lambda1} is provided.}

  \item{min.ratio}{minimal value of infinity-part of the
  penalty that will be tried, as a fraction of the maximal
  \code{lambda1} value. A too small value might lead to
  unstability at the end of the solution path corresponding
  to small \code{lambda1}.  The default value tries to
  avoid this, adapting to the '\eqn{n<p}{n<p}' context.
  Ignored if \code{lambda1} is provided.}

  \item{max.feat}{integer; limits the number of features
  ever to enter the model: in our implementation of the
  bounded regression, it corresponds to the variables which
  have left the boundary along the path.  The algorithm
  stops if this number is exceeded and \code{lambda1} is
  cut at the corresponding level. Default is
  \code{min(nrow(x),ncol(x))} for small \code{lambda2}
  (<0.01) and \code{min(4*nrow(x),ncol(x))} otherwise. Use
  with care, as it considerably changes the computation
  time.}

  \item{control}{list of argument controlling low level
  options of the algorithm --use with care and at your own
  risk-- : \itemize{%

  \item{\code{verbose}: }{integer; activate verbose mode
  --this one is not too much risky!-- set to \code{0} for
  no output; \code{1} for warnings only, and \code{2} for
  tracing the whole progression. Default is \code{1}.
  Automatically set to \code{0} when the method is embedded
  within cross-validation or stability selection.}

  \item{\code{timer}: }{logical; use to record the timing
  of the algorithm. Default is \code{FALSE}.}

  \item{\code{max.iter}: }{the maximal number of iteration
  used to solve the problem for a given value of
  \code{lambda1}. Default is 500.}

  \item{\code{method}: }{a string for the underlying solver
  used. Either \code{"quadra"} or \code{"fista"} are
  available for bounded regression. Default is
  \code{"quadra"}.}

  \item{\code{threshold}: }{a threshold for convergence.
  The algorithm stops when the optimality conditions are
  fulfill up to this threshold. Default is \code{1e-7} for
  \code{"quadra"} and \code{1e-2} for \code{"fista"}.}

  \item{\code{bulletproof}: }{logical; indicates if the
  bulletproof mode should be used while running the
  \code{"quadra"} method. Default is \code{TRUE}. See
  details below.}}}

  \item{checkargs}{logical; should arguments be checked to
  (hopefully) avoid internal crashes? Default is
  \code{TRUE}. Automatically set to \code{FALSE} when calls
  are made from cross-validation or stability selection
  procedures.}
}
\value{
  an object with class \code{quadrupen}, see the
  documentation page \code{\linkS4class{quadrupen}} for
  details.
}
\description{
  Adjust a linear model penalized by a mixture of a
  (possibly weighted) \eqn{\ell_\infty}{l-infinity}-norm
  (bounding the magnitude of the parameters) and a
  (possibly structured) \eqn{\ell_2}{l2}-norm (ridge-like).
  The solution path is computed at a grid of values for the
  infinity-penalty, fixing the amount of \eqn{\ell_2}{l2}
  regularization. See details for the criterion optimized.
}
\note{
  The optimized criterion is \if{latex}{\deqn{%
  \hat{\beta}_{\lambda_1,\lambda_1} = \arg \min_{\beta}
  \frac{1}{2} (y - X \beta)^T (y - X \beta) + \lambda_1 \|D
  \beta \|_{\infty} + \frac{\lambda_2}{2} \beta^T S \beta,
  }} \if{html}{\out{ <center> &beta;<sup>hat</sup>
  <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
  argmin<sub>&beta;</sub> 1/2 RSS(&beta) +
  &lambda;<sub>1</sub> &#124; D &beta;
  &#124;<sub>&infin;</sub> + &lambda;/2 <sub>2</sub>
  &beta;<sup>T</sup> S &beta;, </center> }}
  \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta
  1/2 RSS(beta) + lambda1 max|D beta| + lambda2 beta' S
  beta,}} where \eqn{D}{D} is a diagonal matrix, whose
  diagonal terms are provided as a vector by the
  \code{penscale} argument. The \eqn{\ell_2}{l2}
  structuring matrix \eqn{S}{S} is provided via the
  \code{struct} argument, a positive semidefinite matrix
  (possibly of class \code{Matrix}).

  Note that the quadratic algorithm for the bounded
  regression may become unstable along the path because of
  singularity of the underlying problem, e.g. when there
  are too much correlation or when the size of the problem
  is close to or smaller than the sample size. In such
  cases, it might be a good idea to switch to the proximal
  solver, slower yet more robust. This is the strategy
  adopted by the \code{'bulletproof'} mode, that will send
  a warning while switching the method to \code{'fista'}
  and keep on optimizing on the remainder of the path. When
  \code{bulletproof} is set to \code{FALSE}, the algorithm
  stops at an early stage of the path of solutions. Hence,
  users should be careful when manipulating the resulting
  \code{'quadrupen'} object, as it will not have the size
  expected regarding the dimension of the \code{lambda1}
  argument.

  Singularity of the system can also be avoided with a
  larger \eqn{\ell_2}{l2}-regularization, via
  \code{lambda2}, or a "not-too-small"
  \eqn{\ell_\infty}{l-infinity} regularization, via a
  larger \code{'min.ratio'} argument.
}
\examples{
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) \%*\% chol(Sigma))
y <- 10 + x \%*\% beta + rnorm(n,0,10)

## Infinity norm without/with an additional l2 regularization term
## and with structuring prior
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- "relevant"
plot(bounded.reg(x,y,lambda2=0) , label=labels) ## a mess
plot(bounded.reg(x,y,lambda2=10), label=labels) ## good guys are at the boundaries
plot(bounded.reg(x,y,lambda2=10,struct=solve(Sigma)), label=labels) ## even better
}
\seealso{
  See also \code{\linkS4class{quadrupen}},
  \code{\link{plot.quadrupen}} and \code{\link{crossval}}.
}
\keyword{models,}
\keyword{regression}

