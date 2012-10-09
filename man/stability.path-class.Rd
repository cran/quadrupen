\docType{class}
\name{stability.path-class}
\alias{print,stability.path-method}
\alias{show,stability.path-method}
\alias{stability.path-class}
\title{Class "stability.path"}
\description{
  Class of object returned by the \code{stability}
  function, with methods \code{print}, \code{show} and
  \code{plot}.
}
\section{Slots}{
  \describe{

  \item{\code{probabilities}:}{a \code{Matrix} object
  containing the estimated selection of probabilities along
  the l1-penalty levels.} \item{\code{penalty}:}{Object of
  class \code{"character"} indicating the penalizer used.}
  \item{\code{naive}:}{logical indicating whether rescaling
  of the coefficients has been performed regarding the
  l2-penalty.} \item{\code{lambda1}:}{a vector with the
  l1-penalty levels.} \item{\code{lambda2}:}{a scalar with
  the l2-penalty level.} }
}
\seealso{
  See also \code{\link{plot.stability.path}}, and
  \code{\link{stability}}.
}
\keyword{class}

