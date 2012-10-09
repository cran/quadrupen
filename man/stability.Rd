\name{stability}
\alias{stability}
\title{Stability selection for a quadrupen fit.}
\usage{
  stability(x, y, method = "elastic.net", subsamples = 100,
    sample.size = floor(n/2), randomize = TRUE,
    weakness = 0.5, verbose = TRUE,
    folds = replicate(subsamples, sample(1:nrow(x), sample.size), simplify = FALSE),
    ...)
}
\arguments{
  \item{x}{matrix of features, possibly sparsely encoded
  (experimental). Do NOT include intercept.}

  \item{y}{response vector.}

  \item{method}{a string for the fitting procedure used for
  stability selection. Only \code{"elastic.net"} is
  available for the moment}

  \item{subsamples}{integer indicating the number of
  subsamplings used to estimate the selection
  probabilities. Default is 100.}

  \item{sample.size}{integer indicating the size of each
  subsamples. Default is \code{floor(n/2)}.}

  \item{randomize}{Should a randomized version of the
  fitting procedure by used? Default is \code{TRUE}. See
  details below.}

  \item{weakness}{Coefficient used for randomizing. Default
  is \code{0.5}. Ignored when \code{randomized=TRUE}. See
  details below.}

  \item{folds}{list with \code{subsamples} entries with
  vectors describing the folds to use for the stability
  procedure. By default, the folds are randomly sampled
  with the specified \code{subsamples} argument.}

  \item{verbose}{logical; indicates if the progression
  should be displayed. Default is \code{TRUE}.}

  \item{...}{additional parameters to overwrite the
  defaults of the \code{'method'} fitting procedure. See
  the corresponding documentation (e.g.,
  \code{\link{elastic.net}})}
}
\value{
  An object of class "stability-path".
}
\description{
  Compute the stability path of a (possibly randomized)
  fitting procedure by subsampling as introduced by
  Meinshausen and Buhlmann (2010).
}
\note{
  When \code{randomized = TRUE}, the \code{penscale}
  argument that weights the l1-penalty of the fitting
  procedure is perturbed (divided) for each subsample by a
  random variable uniformly distributed on
  \if{latex}{\eqn{[\alpha,1]}}\if{html}{[&#945;,1]}\if{text}{\eqn{[alpha,1]}},
  alpha being the weakness parameter.
}
\examples{
\dontrun{
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

## Build a vector of label for true nonzeros
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- c("relevant")
labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))

## Call to stability selection function, 200 subsampling
stabout <- stability(x,y, subsamples=200, lambda2=1, min.ratio=1e-2)
## Build the plot an recover the selected variable for a given cutoff
## and per-family error rate
stabpath <- plot(stab, labels=labels, cutoff=0.75, PFER=1)

cat("\\nFalse positives for the randomized Elastic-net with stability selection: ",
     sum(labels[stabpath$selected] != "relevant"))
cat("\\nDONE.\\n")
}
}
\references{
  N. Meinshausen and P. Buhlmann (2010). Stability
  Selection, JRSS(B).
}
\seealso{
  \code{\linkS4class{stability.path}} and
  \code{\link{plot.stability.path}}.
}
\keyword{models,}
\keyword{regression}

