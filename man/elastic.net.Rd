\name{elastic.net}
\alias{elastic.net}
\title{Fit a linear model with elastic-net regularization}
\usage{
  elastic.net(x, y, lambda1 = NULL, lambda2 = 0.01,
    penscale = rep(1, p), struct = NULL, intercept = TRUE,
    normalize = TRUE, naive = FALSE,
    nlambda1 = ifelse(is.null(lambda1), 100, length(lambda1)),
    min.ratio = ifelse(n <= p, 0.01, 1e-04),
    max.feat = ifelse(lambda2 < 0.01, min(n, p), min(4 * n, p)),
    beta0 = NULL, control = list(), checkargs = TRUE)
}
\arguments{
  \item{x}{matrix of features, possibly sparsely encoded
  (experimental). Do NOT include intercept. When normalized
  os \code{TRUE}, coefficients will then be rescaled to the
  original scale.}

  \item{y}{response vector.}

  \item{lambda1}{sequence of decreasing
  \eqn{\ell_1}{l1}-penalty levels. If \code{NULL} (the
  default), a vector is generated with \code{nlambda1}
  entries, starting from a guessed level \code{lambda1.max}
  where only the intercept is included, then shrunken to
  \code{min.ratio*lambda1.max}.}

  \item{lambda2}{real scalar; tunes the \eqn{\ell_2}{l2}
  penalty in the Elastic-net. Default is 0.01. Set to 0 to
  recover the Lasso.}

  \item{penscale}{vector with real positive values that
  weight the \eqn{\ell_1}{l1}-penalty of each feature.
  Default set all weights to 1.}

  \item{struct}{matrix structuring the coefficients,
  possibly sparsely encoded. Must be at least positive
  semidefinite (this is checked internally if the
  \code{checkarg} argument is \code{TRUE}). If \code{NULL}
  (the default), the identity matrix is used. See details
  below.}

  \item{intercept}{logical; indicates if an intercept
  should be included in the model. Default is \code{TRUE}.}

  \item{normalize}{logical; indicates if variables should
  be normalized to have unit L2 norm before fitting.
  Default is \code{TRUE}.}

  \item{naive}{logical; Compute either 'naive' of classic
  elastic-net as defined in Zou and Hastie (2006): the
  vector of parameters is rescaled by a coefficient
  \code{(1+lambda2)} when \code{naive} equals \code{FALSE}.
  No rescaling otherwise. Default is \code{FALSE}.}

  \item{nlambda1}{integer that indicates the number of
  values to put in the \code{lambda1} vector.  Ignored if
  \code{lambda1} is provided.}

  \item{min.ratio}{minimal value of \eqn{\ell_1}{l1}-part
  of the penalty that will be tried, as a fraction of the
  maximal \code{lambda1} value. A too small value might
  lead to unstability at the end of the solution path
  corresponding to small \code{lambda1} combined with
  \eqn{\lambda_2=0}{lambda2=0}.  The default value tries to
  avoid this, adapting to the '\eqn{n<p}{n<p}' context.
  Ignored if \code{lambda1} is provided.}

  \item{max.feat}{integer; limits the number of features
  ever to enter the model; i.e., non-zero coefficients for
  the Elastic-net: the algorithm stops if this number is
  exceeded and \code{lambda1} is cut at the corresponding
  level. Default is \code{min(nrow(x),ncol(x))} for small
  \code{lambda2} (<0.01) and \code{min(4*nrow(x),ncol(x))}
  otherwise. Use with care, as it considerably changes the
  computation time.}

  \item{beta0}{a starting point for the vector of
  parameter. When \code{NULL} (the default), will be
  initialized at zero. May save time in some situation.}

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
  used to solve the problem for a given value of lambda1.
  Default is 500.}

  \item{\code{method}: }{a string for the underlying solver
  used. Either \code{"quadra"}, \code{"pathwise"} or
  \code{"fista"}. Default is \code{"quadra"}.}

  \item{\code{threshold}: }{a threshold for convergence.
  The algorithm stops when the optimality conditions are
  fulfill up to this threshold. Default is \code{1e-7} for
  \code{"quadra"} and \code{1e-2} for the first order
  methods.}

  \item{\code{monitor}: }{indicates if a monitoring of the
  convergence should be recorded, by computing a lower
  bound between the current solution and the optimum: when
  \code{'0'} (the default), no monitoring is provided; when
  \code{'1'}, the bound derived in Grandvalet et al. is
  computed; when \code{'>1'}, the Fenchel duality gap is
  computed along the algorithm.} }}

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
  Adjust a linear model with elastic-net regularization,
  mixing a (possibly weighted) \eqn{\ell_1}{l1}-norm
  (LASSO) and a (possibly structured) \eqn{\ell_2}{l2}-norm
  (ridge-like). The solution path is computed at a grid of
  values for the \eqn{\ell_1}{l1}-penalty, fixing the
  amount of \eqn{\ell_2}{l2} regularization. See details
  for the criterion optimized.
}
\note{
  The optimized criterion is the following:
  \if{latex}{\deqn{% \hat{\beta}_{\lambda_1,\lambda_1} =
  \arg \min_{\beta} \frac{1}{2} (y - X \beta)^T (y - X
  \beta) + \lambda_1 \|D \beta \|_{1} + \frac{\lambda_2}{2}
  \beta^T S \beta, }} \if{html}{\out{ <center>
  &beta;<sup>hat</sup>
  <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
  argmin<sub>&beta;</sub> 1/2 RSS(&beta) +
  &lambda;<sub>1</sub> &#124; D &beta; &#124;<sub>1</sub> +
  &lambda;/2 <sub>2</sub> &beta;<sup>T</sup> S &beta;,
  </center> }} \if{text}{\deqn{beta.hat(lambda1, lambda2) =
  argmin_beta 1/2 RSS(beta) + lambda1 |D beta|1 + lambda2
  beta' S beta,}} where \eqn{D}{D} is a diagonal matrix,
  whose diagonal terms are provided as a vector by the
  \code{penscale} argument. The \eqn{\ell_2}{l2}
  structuring matrix \eqn{S}{S} is provided via the
  \code{struct} argument, a positive semidefinite matrix
  (possibly of class \code{Matrix}).
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

## Structured Elastic.net
p <- ncol(x)
C <- bandSparse(p,k=0:1,diagonals=list(rep(1,p),rep(-1,p-1)))

## This gives a great advantage to the elastic-net
## for support recovery
beta.lasso <- slot(crossval(x,y,lambda2=0 , mc.cores=2) , "beta.min")
beta.enet  <- slot(crossval(x,y,lambda2=10, mc.cores=2), "beta.min")
beta.struc <- slot(crossval(x,y,lambda2=10, mc.cores=2, struct=t(C) \%*\% C), "beta.min")

cat("\\nFalse positives for the Lasso:", sum(sign(beta) != sign(beta.lasso)))
cat("\\nFalse positives for the Elastic-net:", sum(sign(beta) != sign(beta.enet)))
cat("\\nDONE.\\n")

labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- "relevant"
## Comparing the solution path pf the LASSO, the Elastic-net and the
## Structured Elastic-net
plot(elastic.net(x,y,lambda2=0), label=labels) ## a mess
plot(elastic.net(x,y,lambda2=10), label=labels) ## a lot better
plot(elastic.net(x,y,lambda2=10,struct=solve(Sigma)), label=labels) ## even better
}
\seealso{
  See also \code{\linkS4class{quadrupen}},
  \code{\link{plot.quadrupen}} and \code{\link{crossval}}.
}
\keyword{models,}
\keyword{regression}

