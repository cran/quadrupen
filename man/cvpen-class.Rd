\docType{class}
\name{cvpen-class}
\alias{cvpen-class}
\title{Class "cvpen"}
\description{
  Class of object returned by a cross-validation performed
  through the \code{crossval} method.
}
\section{Slots}{
  \describe{ \item{\code{lambda1}:}{vector of
  \eqn{\lambda_1}{lambda1} (\eqn{\ell_1}{l1} or
  \eqn{\ell_\infty}{infinity} penalty levels) for which
  each cross-validation has been performed.}
  \item{\code{lambda2}:}{vector (or scalar) of
  \eqn{\ell_2}{l2}-penalty levels for which each
  cross-validation has been performed.}
  \item{\code{lambda1.min}:}{level of
  \eqn{\lambda_1}{lambda1} that minimizes the error
  estimated by cross-validation.}
  \item{\code{lambda1.1se}:}{largest level of
  \eqn{\lambda_1}{lambda1} such as the cross-validated
  error is within 1 standard error of the minimum.}
  \item{\code{lambda2.min}:}{level of
  \eqn{\lambda_2}{lambda2} that minimizes the error
  estimated by cross-validation.} \item{\code{cv.error}:}{a
  data frame containing the mean cross-validated error and
  its associated standard error for each values of
  \code{lambda1} and \code{lambda2}.}
  \item{\code{folds}:}{list of \code{K} vectors indicating
  the folds used for cross-validation.}
  \item{\code{beta.min}:}{the vector of parameters obtained
  by fitting the problem on the full data set \code{x} and
  \code{y} with \code{lambda1.min} and \code{lambda2.min}
  penalties.} \item{\code{beta.1se}:}{the vector of
  parameters obtained by fitting the problem on the full
  data set \code{x} and \code{y} with \code{lambda1.1se}
  and \code{lambda2.min} penalties. }

  }

  The specific \code{\link{plot.cvpen}} method is
  documented.
}
\seealso{
  See also \code{\link{plot.cvpen}} and
  \code{\link{crossval}}.
}
\keyword{class}

