\docType{methods}
\name{plot,stability.path-method}
\alias{plot,stability.path-method}
\alias{plot.stability.path}
\title{Plot method for \code{stability.path}.}
\arguments{
  \item{x}{output of a \code{crossval} run (must be of
  class \code{stability.path}).}

  \item{y}{used for S4 compatibility.}

  \item{main}{main title.}

  \item{xvar}{variable to plot on the X-axis: either
  \code{"ave.sel"} (average number of selected features),
  \code{"lambda"} (l1-penalty level) or \code{"fraction"}
  (fraction of l1-penalty level applied). Default is set to
  \code{"ave.sel"}.}

  \item{log.scale}{logical; indicates if a log-scale should
  be used when \code{xvar="lambda"}.}

  \item{plot}{logical; indicates if the graph should be
  plotted. Default is \code{TRUE}.}

  \item{annot}{logical; should annotation be made on the
  graph regarding controlled PFER. Default is \code{TRUE}.}

  \item{cutoff}{value of the cutoff probability. See
  details.}

  \item{PFER}{value of the per-family error rate to
  control. See details.}

  \item{labels}{an optional vector of labels for each
  variable in the path (e.g., 'relevant'/'irrelevant').}

  \item{...}{used for S4 compatibility.}
}
\value{
  a list with a ggplot2 object which can be plotted via the
  \code{print} method, and a vector of selected variables
  corresponding to the given cutoff and PFER argument.
}
\description{
  Produce a plot of the stability path obtained by
  stability selection. Display the average number of
  selected variables to controlled the per-family error
  rate at a given cutoff probability.
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
stabpath <- plot(stab, xvar="lambda", labels=labels, cutoff=0.75, PFER=2)
stabpath <- plot(stab, xvar="fraction", annot=FALSE, labels=labels, cutoff=0.75, PFER=1)
}
}

