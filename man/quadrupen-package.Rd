\docType{package}
\name{quadrupen-package}
\alias{quadrupen}
\alias{quadrupen-package}
\title{Sparsity by Worst-Case Quadratic Penalties}
\description{
  This package is designed to fit accurately several
  popular penalized linear regression models using the
  algorithm proposed in Grandvalet et al. (submitted) by
  solving quadratic problems with increasing size.
}
\section{Features}{
  For now, an unique fitting function \code{elastic.net} is
  available, which, depending on how the parameters are
  tuned, can solve \itemize{%

  \item the LASSO of Tibshirani (1996), \item the
  adaptive-LASSO of Zou (2006), \item the Elastic-net of
  Zou and Hastie (2006), \item the Structured Elastic-net
  of Slawsky et al. (2010). }

  The S4 objects produced by the fitting procedures own the
  classical methods for linear model in \code{R}, as well
  as methods for plotting, cross-validating, and for the
  stability selection (Meinshausen and Buhlmann, 2010). See
  the available demo through \code{demo(quad_enet)} and
  examples.
}

\section{Algorithm}{
  The general strategy of the algorithm relies on
  maintaining an active set of variables, starting from a
  vector of zeros. The underlying optimization problem is
  solved only on the activated variables, thus handling
  with small smooth problems with increasing size. The
  active set algorithm stops when the optimality conditions
  are met.

  Three strategies are implemented to successively solve
  the underlying optimization problem \itemize{%

  \item the direct resolution of quadratic sub problems, as
  described in Grandvalet et al. (submitted), \item a
  coordinate descent approach, \item a proximal (FISTA)
  approach. }

  The default settings uses the quadratic approach that
  gives its name to the package. It has been optimized to
  be the method of choice for small and medium scale
  problems, and produce very accurate solutions.
  Nevertheless, the first order methods (coordinate descent
  and FISTA) can be interesting in situations at the end of
  the solution path where the problem is close to singular,
  in which case the Cholesky decomposition used in the
  quadratic solver can be computationally unstable.
}

\section{Technical remarks}{
  Most of the algebra are performed in C++, relying on the
  \pkg{RcppArmadillo} package. We also provide a (double)
  cross-validation procedure which uses the multi-core
  capability of the computer, through the \pkg{parallel}
  package. This feature is only available for Windows user,
  though. Finally, note that the plot methods enjoy some
  (still very few) of the capabilities of \pkg{ggplot2}
  package.

  We hope to enrich the package with other popular fitting
  procedures and provide other tools, particularly towards
  bootstrapping and model selection purpose. Sparse matrix
  encoding is partially supported for the moment, and will
  hopefully be fully available in the future.
}
\author{
  \email{julien.chiquet@genopole.cnrs.com}
}
\references{
  \itemize{%

  \item Yves Grandvalet, Christophe Ambroise, Julien
  Chiquet Sparsity by Worst-case Quadratic Penalties, arXiv
  preprint. \item Nicolas Meinshausen and Peter Buhlmann.
  Stability Selection, JRSS(B), 2010. \item Martin Slawski,
  Wolfgang zu Castell, and Gerhard Tutz. Feature selection
  guided by structural information, AOAS, 2010. \item
  Robert Tibshirani. Regression Shrinkage and Selection Via
  the Lasso, JRSS(B), 1996. \item Hui Zou. The Adaptive
  Lasso and Its Oracle Properties, JASA, 2006. \item Hui
  Zou and Trevor Hastie. Regularization and variable
  selection via the elastic net, JRSS(B), 2006. }
}

