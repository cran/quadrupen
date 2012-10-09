\docType{methods}
\name{plot,quadrupen-method}
\alias{plot,quadrupen-method}
\alias{plot.quadrupen}
\title{Plot method for a quadrupen object}
\arguments{
  \item{x}{output of a fitting procedure of the
  \pkg{quadrupen} package (only \code{\link{elastic.net}}
  for the moment). Must be of class \code{quadrupen}.}

  \item{y}{used for S4 compatibility.}

  \item{xvar}{variable to plot on the X-axis: either
  \code{"lambda"} (l1-penalty level) or \code{"fraction"}
  (l1-norm of thje coefficients). Default is set to
  \code{"lambda"}.}

  \item{main}{the main title. Default is set to the model
  name followed by what is on the Y-axis.}

  \item{log.scale}{logical; indicates if a log-scale should
  be used when \code{xvar="lambda"}.}

  \item{standardize}{logical; standardize the coefficient
  before plotting. Default is \code{TRUE}.}

  \item{reverse}{logical; should the X-axis be reversed
  when \code{xvar="lambda"}? Default is \code{FALSE}.}

  \item{label}{vector indicating the names associated to
  the plotted variables. When specified, a legend is drawn
  in order to identify each variable. Only relevant when
  the number of predictor is small. Remind that the
  intercept does not count. Default is \code{NULL}.}

  \item{plot}{logical; indicates if the graph should be
  plotted on call. Default is \code{TRUE}.}
}
\value{
  a ggplot2 object which can be plotted via the
  \code{print} method.
}
\description{
  Produce a plot of the solution path of a \code{quadrupen}
  fit.
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
## This gives a great advantage to the elastic-net
## for support recovery; have a look at image(Sigma)
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) \%*\% chol(Sigma))
y <- 10 + x \%*\% beta + rnorm(n,0,10)
lasso <- elastic.net(x,y, lambda2=0)
enet  <- elastic.net(x,y, lambda2=10)

## Plot the Lasso path
plot(lasso, main="Lasso solution path")
## Plot the Elastic-net path
plot(enet, main = "Elastic-net solution path")
## Plot the Elastic-net path (fraction on X-axis, unstandardized coefficient)
plot(enet, standardize=FALSE, xvar="fraction")
}
}
\seealso{
  \code{\linkS4class{quadrupen}}.
}

