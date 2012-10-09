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
  \pkg{quadrupen} package, e.g.; \code{elastic.net}.
}
\section{Slots}{
  \describe{

  \item{\code{coefficients}:}{Matrix (class
  \code{"Matrix"}) of coefficients with respect to the
  original input (even when normalized). The number of rows
  of \code{coefficients} corresponds the length of
  \code{lambda1}.}

  \item{\code{intercept}:}{logical; indicates if an
  intercept has been included to the model.}

  \item{\code{normalize}:}{logical; indicates if the matrix
  of features has been normalized for the fitting.}

  \item{\code{mu}:}{A vector (class \code{"numeric"})
  containing the successive value of the (unpenalized)
  intercept. Equals to zero if \code{intercept} has been
  set to \code{FALSE}.}

  \item{\code{meanx}:}{Vector (class \code{"numeric"})
  containing the mean of the design ' matrix. Vector of
  zero if \code{intercept} has been set to \code{FALSE}.}

  \item{\code{normx}:}{Vector (class \code{"numeric"})
  containing the square root of the sum of squares of each
  columns of the design matrix. Vector of one if normalized
  has been set to \code{FALSE}.}

  \item{\code{penscale}:}{Vector \code{"numeric"} with real
  positive values that weight the l1-penalty of each
  feature.}

  \item{\code{penalty}:}{Object of class \code{"character"}
  indicating the penalizer used (only "elastic-net"
  available for the moment).}

  \item{\code{naive}:}{logical; indicates whether the
  'naive' or the 'classic' elastic-net has been fitted.}

  \item{\code{lambda1}:}{Vector (class \code{"numeric"}) of
  l1-penalty levels for which the model has eventually been
  fitted.}

  \item{\code{lambda2}:}{Scalar (class \code{"numeric"})
  for the amount of ridge-like penalty.}

  \item{\code{struct}:}{\code{"Matrix"} structuring the
  coefficients, possibly sparsely encoded (MUST be positive
  definite).}

  \item{\code{control}:}{Object of class \code{"list"}.}

  \item{\code{monitoring}:}{List (class \code{"list"})
  which contains various indicator dealing with the
  optimization process.}

  \item{\code{residuals}:}{Vector (class \code{"numeric"})
  of (weighted) residuals.}

  \item{\code{fitted}:}{Vector (class \code{"numeric"}) of
  fitted values.}

  \item{\code{weights}:}{Vector of observation weights
  used.} }
}

\section{Methods}{
  This class comes with the usual \code{predict(object,
  newx, ...)}, \code{fitted(object, ...)},
  \code{residuals(object, ...)}, \code{print(object, ...)},
  \code{show(object)} and \code{deviance(object, ...)}
  generic (undocumented) methods.

  Two specific and documented methods are also available:

  \describe{ \item{plot}{plots a \code{quadrupen} object}
  \item{crossval}{cross-validates a \code{quadrupen}
  object} }
}
\seealso{
  See also \code{\link{plot.quadrupen}}, and
  \code{\link{crossval}}.
}
\keyword{class}

