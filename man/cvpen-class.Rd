\docType{class}
\name{cvpen-class}
\alias{cvpen-class}
\title{Class "cvpen"}
\description{
  Class of object returned by a cross-validation performed
  through the \code{crossval} method.
}
\section{Slots}{
  \describe{ \item{\code{lambda1}:}{vector of l1-penalty
  levels for which each cross-validation has been
  performed.} \item{\code{lambda2}:}{vector (or scalar) of
  l2-penalty levels for which each cross-validation has
  been performed.} \item{\code{lambda1.min}:}{level of
  l1-penalty that minimizes the error estimated by
  cross-validation.} \item{\code{lambda1.1se}:}{largest
  level of l1-penalty such has the cross-validated error is
  within 1 standard error of the minimum.}
  \item{\code{lambda2.min}:}{level of l2-penalty that
  minimizes the error estimated by cross-validation.}
  \item{\code{cv.error}:}{a data frame containing the mean
  cross-validated error and its associated standard error
  for each values of \code{lamdba1} and \code{lamdba2}.}
  \item{\code{folds}:}{list of \code{K} vectors indicating
  the folds used for cross-validation.}
  \item{\code{beta.min}:}{the vector of parameters obtained
  by fitting the problem on the full data set \code{x} and
  \code{y} with l1-penalty \code{lambda1.min} and
  l2-penalty \code{lambda2.min}.}
  \item{\code{beta.1se}:}{the vector of parameters obtained
  by fitting the problem on the full data set \code{x} and
  \code{y} with l1-penalty \code{lambda1.1se} and
  l2-penalty \code{lambda2.min}. }

  }

  The specific \code{\link{plot.cvpen}} method is
  documented.
}
\seealso{
  See also \code{\link{plot.cvpen}} and
  \code{\link{crossval}}.
}
\keyword{class}

