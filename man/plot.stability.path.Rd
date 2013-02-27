\docType{methods}
\name{plot,stability.path-method}
\alias{plot.stability.path}
\alias{plot,stability.path-method}
\title{Plot method for \code{stability.path}.}
\arguments{
  \item{x}{output of a \code{stability} run (must be of
  class \code{stability.path}).}

  \item{y}{used for S4 compatibility.}

  \item{main}{main title. If none given, a somewhat
  appropriate title is automatically generated.}

  \item{xvar}{variable to plot on the X-axis: either
  \code{"lambda"} (first penalty level) or
  \code{"fraction"} (fraction of the penalty level applied
  tune by \eqn{\lambda_1}{lambda1}). Default is
  \code{"lambda"}.}

  \item{log.scale}{logical; indicates if a log-scale should
  be used when \code{xvar="lambda"}. Default is
  \code{TRUE}.}

  \item{plot}{logical; indicates if the graph should be
  plotted. Default is \code{TRUE}. If \code{FALSE}, only
  the \pkg{ggplot2} object is sent back.}

  \item{sel.mode}{a character string, either \code{'rank'}
  or \code{'PFER'}. In the first case, the selection is
  based on the rank of total probabilties by variables
  along the path: the first \code{nvar} variables are
  selected (see below). In the second case, the PFER
  control is used as described in Meinshausen and
  Buhlmannn's paper. Default is \code{'rank'}.}

  \item{nvar}{number of variables selected (only relevant
  when \code{sel.mode} equals \code{'rank'}. Default is
  \code{floor(n/log(p))}.}

  \item{cutoff}{value of the cutoff probability (only
  relevant when \code{sel.mode} equals \code{'PFER'}).}

  \item{PFER}{value of the per-family error rate to control
  (only relevant when \code{sel.mode} equals
  \code{'PFER'}).}

  \item{labels}{an optional vector of labels for each
  variable in the path (e.g., 'relevant'/'irrelevant'). See
  examples.}

  \item{annot}{logical; should annotation be made on the
  graph regarding controlled PFER (only relevant when
  \code{sel.mode} equals \code{'PFER'})? Default is
  \code{TRUE}.}

  \item{...}{used for S4 compatibility.}
}
\value{
  a list with a \pkg{ggplot2} object which can be plotted
  via the \code{print} method, and a vector of selected
  variables corresponding to method of choice
  (\code{'rank'} or \code{'PFER'}).
}
\description{
  Produce a plot of the stability path obtained by
  stability selection.
}
\examples{
\dontrun{
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
stab <- stability(x,y, subsamples=200, lambda2=1, min.ratio=1e-2)

## Build the plot an recover the selected variable
plot(stab, labels=labels)
plot(stab, xvar="fraction", labels=labels, sel.mode="PFER", cutoff=0.75, PFER=2)
}
}

