\docType{class}
\name{quadrupen-class}
\alias{deviance,quadrupen-method}
\alias{fitted,quadrupen-method}
\alias{predict,quadrupen-method}
\alias{print,quadrupen-method}
\alias{quadrupen-class}
\alias{residuals,quadrupen-method}
\alias{show,quadrupen-method}
\title{Class "quadrupen"}
\description{
  Class of object returned by any fitting function of the
  \pkg{quadrupen} package (\code{elastic.net} or
  \code{bounded.reg}).
}
\section{Slots}{
  \describe{

  \item{\code{coefficients}:}{Matrix (class
  \code{"dgCMatrix"}) of coefficients with respect to the
  original input. The number of rows corresponds the length
  of \code{lambda1}.}

  \item{\code{active.set}:}{Matrix (class
  \code{"dgCMatrix"}, generally sparse) indicating the
  'active' variables, in the sense that they activate the
  constraints. For the \code{\link{elastic.net}}, it
  corresponds to the nonzero variables; for the
  \code{\link{bounded.reg}} function, it is the set of
  variables reaching the boudary along the path of
  solutions.}

  \item{\code{intercept}:}{logical; indicates if an
  intercept has been included to the model.}

  \item{\code{mu}:}{A vector (class \code{"numeric"})
  containing the successive values of the (unpenalized)
  intercept. Equals to zero if \code{intercept} has been
  set to \code{FALSE}.}

  \item{\code{meanx}:}{Vector (class \code{"numeric"})
  containing the column means of the predictor matrix.}

  \item{\code{normx}:}{Vector (class \code{"numeric"})
  containing the square root of the sum of squares of each
  column of the design matrix.}

  \item{\code{penscale}:}{Vector \code{"numeric"} with real
  positive values that have been used to weight the penalty
  tuned by \eqn{\lambda_1}{lambda1}.}

  \item{\code{penalty}:}{Object of class \code{"character"}
  indicating the method used (\code{"elastic-net"} or
  \code{"bounded regression"}).}

  \item{\code{naive}:}{logical; was the \code{naive} mode
  on?}

  \item{\code{lambda1}:}{Vector (class \code{"numeric"}) of
  penalty levels (either \eqn{\ell_1}{l1} or
  \eqn{\ell_\infty}{l-infinity}) for which the model has
  eventually been fitted.}

  \item{\code{lambda2}:}{Scalar (class \code{"numeric"})
  for the amount of \eqn{\ell_2}{l2} (ridge-like) penalty.}

  \item{\code{struct}:}{Object of class \code{"Matrix"}
  used to structure the coefficients in the
  \eqn{\ell_2}{l2} penalty.}

  \item{\code{control}:}{Object of class \code{"list"} with
  low level options used for optimization.}

  \item{\code{monitoring}:}{List (class \code{"list"})
  which contains various indicators dealing with the
  optimization process.}

  \item{\code{residuals}:}{Matrix of residuals, each column
  corresponding to a value of \code{lambda1}.}

  \item{\code{r.squared}:}{Vector (class \code{"numeric"})
  given the coefficient of determination as a function of
  lambda1.}

  \item{\code{fitted}:}{Matrix of fitted values, each
  column corresponding to a value of \code{lambda1}.} }
}

\section{Methods}{
  This class comes with the usual \code{predict(object,
  newx, ...)}, \code{fitted(object, ...)},
  \code{residuals(object, ...)}, \code{print(object, ...)},
  \code{show(object)} and \code{deviance(object, ...)}
  generic (undocumented) methods.

  A specific plotting method is available and documented
  (\code{\link{plot.quadrupen}}).
}
\seealso{
  See also \code{\link{plot.quadrupen}}.
}
\keyword{class}

