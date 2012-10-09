##' Fit a linear model with elastic-net regularization
##' 
##' Adjust a linear model with elastic-net regularization, mixing a
##' (possibly weighted) l1-penalty (LASSO) and a (possibly structured)
##' l2-penalty (ridge-like). The solution path is computed at a grid
##' of values for the l1-penalty, fixing the value of the
##' l2-penalty. See details for the criteria optimized.
##' 
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept.
##' 
##' @param y response vector.
##' 
##' @param lambda1 sequence of decreasing l1-penalty levels. If
##' \code{NULL} (the default), a vector is generated with
##' \code{nlambda1} entries, starting from a guessed level
##' \code{lambda1.max} where only the intercept is included, then
##' shrunken to \code{min.ratio*lambda1.max}.
##' 
##' @param lambda2 real scalar; tunes the ridges penalty in the
##' Elastic-net. Default is 0.01. Set to 0 to recover the Lasso.
##' 
##' @param penscale vector with real positive values that weight the
##' l1-penalty of each features. Default set all weights to 1.
##' 
##' @param struct matrix structuring the coefficients, possibly
##' sparsely encoded (MUST be positive definite). Default is the
##' identity matrix. See details below.
##' 
##' @param weights observation weights for weighted
##' least-squares. Default is 1 for each observation.
##' 
##' @param intercept logical; indicates if an intercept should be
##' included in the model. Default is \code{TRUE}.
##' 
##' @param normalize logical; indicates if variables should be
##' normalized to have unit L2 norm before fitting.  Coefficients will
##' then be rescaled to the original scale. Default is \code{TRUE}.
##'
##' @param naive logical; Compute either 'naive' of classic
##' elastic-net as defined in Zou and Hastie (2006). Default is
##' FALSE. See details below.
##' 
##' @param nlambda1 integer that indicates the number of values to put
##' in the \code{lambda1} vector.  Ignored if \code{lambda1} is
##' provided.
##' 
##' @param min.ratio minimal value of l1-penalty that will be tried,
##' as a fraction of the maximal \code{lambda1} value. A too small
##' value of \code{lambda1} might lead to unstable path of solution
##' when n<p.  The default value tries to avoid this, adapting to the
##' 'n<p' context. Ignored if \code{lambda1} is provided.
##' 
##' @param max.feat integer; limits the number of features ever to
##' enter the model: the algorithm stops if this number is exceeded
##' and \code{lambda1} is cut at the corresponding level. Default is
##' \code{min(nrow(x),ncol(x))} for small \code{lambda2} (<0.01) and
##' \code{min(3*nrow(x),ncol(x))} otherwise. Use with care, as it
##' considerably changes the computation time.
##' 
##' @param control list of argument controlling low level options of
##' the algorithm --use with care and at your own risk-- :
##' \itemize{%
##' 
##' \item{\code{verbose}:}{logical; activate verbose mode. Default is \code{FALSE}.}
##' \item{\code{timer}:}{logical; use to record the timing of the algorithm. Default is \code{FALSE}.}
##' \item{\code{zero}:}{a practical zero. Default is \code{.Machine$double.eps}}
##' \item{\code{max.iter}:}{the maximal number of iteration used to solve the problem for a given value of lambda1. Default is 500.}
##' \item{\code{method}:}{a string for the underlying solver used. Either \code{"quadra"}, \code{"pathwise"} or \code{"fista"}. Default is \code{"quadra"}}
##'  \item{\code{threshold}:}{a threshold for convergence. The algorithm stops
##' when the optimality conditions are fulfill up to this
##' threshold. Default is \code{1e-7} for \code{"quadra"} and
##' \code{1e-2} for the first order methods.}
##' \item{\code{monitor}:}{indicates if a monitoring of the convergence should
##' be recorded, by computing a lower bound between the current
##' solution and the optimum: when \code{'0'} (the default), no
##' monitoring is provided; when \code{'1'}, the bound derived in
##' Grandvalet et al. is computed; when \code{'>1'}, the Fenchel
##' duality gap is computed along the algorithm.}}
##'
##' @param checkargs logical; should arguments be checked to
##' (hopefully) avoid internal crashes? Default is \code{TRUE}.
##' 
##' @return a object with class \code{quadrupen}, see the
##' documentation page \code{\linkS4class{quadrupen}} for details.
##' 
##' @note The optimized criterion is the following:
##' \if{latex}{\deqn{% \hat{\beta}_{\lambda_1,\lambda_1} = \arg
##' \min_{\beta} \frac{1}{2} (y - X \beta)^T D (y - X \beta) +
##' \lambda_1 \|w \circ \beta \|_1 + \frac{\lambda_2}{2} \beta^T S
##' \beta, }} \if{html}{\out{ <center> &beta;<sup>hat</sup>
##' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
##' argmin<sub>&beta;</sub> 1/2 RSS(&beta;,weights) +
##' &lambda;<sub>1</sub> &#124; &beta; &#124;<sub>1</sub> + &lambda;/2
##' <sub>2</sub> &beta;<sup>T</sup> S &beta;, </center> }}
##' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
##' RSS(beta,weights) + lambda1 ||w.beta||1 + lambda2 beta' S beta,}}
##' where the vector w is the \code{penscale} argument and S is the
##' \code{struct} argument, a positive definite matrix (possibly of
##' class \code{Matrix}). \if{latex}{\eqn{D} is a diagonal matrix
##' whose diagonal corresponds to the 'weights' vector to perform
##' weighted least-squares.}
##'
##' If \code{naive} is \code{TRUE}, the coefficients are not
##' rescale. When \code{naive} is \code{FALSE} (the default) their are
##' rescale by a scalar \code{1+lambda2}, as suggested in Zou and
##' Hastie, to obtained the so-called Elastic-net estimator.
##' 
##' @seealso See also \code{\linkS4class{quadrupen}},
##' \code{\link{plot.quadrupen}} and \code{\link{crossval}}.
##' @name elastic.net
##' @rdname elastic.net
##' @keywords models, regression
##'
##' @examples 
##' rm(list=ls())
##' library(quadrupen)
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
##' Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
##' diag(Sigma) <- 1
##' n <- 100
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##' 
##' ## This gives a great advantage to the elastic-net
##' ## for support recovery
##' beta.lasso <- slot(crossval(x,y,lambda2=0) , "beta.min")
##' beta.enet  <- slot(crossval(x,y,lambda2=10), "beta.min")
##' 
##' cat("\nFalse positives for the Lasso:", sum(sign(beta) != sign(beta.lasso)))
##' cat("\nFalse positives for the Elastic-net:", sum(sign(beta) != sign(beta.enet)))
##' cat("\nDONE.\n")
##' 
##' @export
elastic.net <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.01,
                        penscale  = rep(1,p),
                        weights   = rep(1,n),
                        struct    = 1*bandSparse(p, k=0),
                        intercept = TRUE,
                        normalize = TRUE,
                        naive     = FALSE,
                        nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                        min.ratio = ifelse(n<p,0.01,5e-3),
                        max.feat  = as.integer(ifelse(lambda2<1e-2,min(n,p),min(3*n,p))),
                        control   = list(),
                        checkargs = TRUE) {
  
  ## ===================================================
  ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
  p <- ncol(x) # problem size
  n <- nrow(x) # sample size
  if (checkargs) {
    if(!inherits(x, c("matrix", "Matrix")))
      stop("x has to be of class 'matrix', 'dgeMatrix' or 'dgCMatrix'.")
    if(any(is.na(x)))
      stop("NA value in x not allowed.")
    if(!is.numeric(y))
      stop("y has to be of type 'numeric'")
    if(n != length(y))
      stop("x and y have not correct dimensions")
    if(length(penscale) != p)
      stop("penscale must have ncol(x) entries")
    if (any(penscale <= 0))
      stop("weights in penscale must be positive")
    if(!inherits(lambda2, "numeric") | length(lambda2) > 1)
      stop("lambda2 must be a scalar.")
    if(lambda2 < 0)
      stop("lambda2 must be a non negative scalar.")
    if (!is.null(lambda1)) {
      if(any(lambda1 <= 0))
        stop("entries inlambda1 must all be postive.")
      if(is.unsorted(rev(lambda1)))
        stop("lambda1 values must be sorted in decreasing order.")
    }
    if(min.ratio < 0)
      stop("min.ratio must be non negative.")
    if (ncol(struct) != p | ncol(struct) != p)
      stop("struct must be a (square) positive definite matrix.")
    if (length(max.feat)>1)
      stop("max.feat must be an integer.")
    if(is.numeric(max.feat) & !is.integer(max.feat))
      max.feat <- as.integer(max.feat)
    if (length(weights) != n)
      stop("weights must have n entries.")
    if (any(weights < 0))
      stop("weights must all be non negative.")
  }
  ## ============================================
  ## RECOVERING LOW LEVEL OPTIONS
  quadra <- TRUE
  if (!is.null(control$method)) {
    if (control$method != "quadra") {
      quadra <- FALSE
    }
  }
  ctrl <- list(verbose   = FALSE, # default control options
               timer     = FALSE,
               zero      = .Machine$double.eps,
               max.iter  = 500,
               method    = "quadra",
               threshold = ifelse(quadra, 1e-7, 1e-2),
               monitor   = 0)
  ctrl[names(control)] <- control # overwrite by user specifications
  if (ctrl$timer) {r.start <- proc.time()}

  ## ======================================================
  ## INTERCEPT AND NORMALIZATION TREATMENT
  input <- standardize(x,y,weights,intercept,normalize,ctrl$zero)  

  ## ======================================================
  ## GENERATE A GRID OF PENALTY IF NONE HAS BEEN PROVIDED
  if (is.null(lambda1))
    lambda1 <- get.lambda1(input$xty,penscale,nlambda1,min.ratio)
  
  ## ======================================================
  ## STRUCTURATING MATRIX
  if (lambda2 > 0) {
    struct <- as(struct, "dgCMatrix")
    S <- list(Si = struct@i, Sp = struct@p, Sx = lambda2*struct@x)
  } else {
    S <- NULL
  }
  
  ## ======================================================
  ## STARTING C++ CALL TO ENET_LS
  if (ctrl$timer) {cpp.start <- proc.time()}
  out <- .Call("elastic_net",
               input$x      ,
               input$xty    ,
               S            ,
               lambda1      ,
               penscale     ,
               input$xbar   ,
               input$ybar   ,
               input$normx  ,
               input$normy  ,
               weights      ,
               naive        ,
               ctrl$thresh  ,
               ctrl$zero    ,
               ctrl$max.iter,
               max.feat     ,
               switch(ctrl$method,
                      quadra   = 0,
                      pathwise = 1,
                      fista    = 2, 0),
               ctrl$verbose,
               inherits(x, "sparseMatrix"),
               ctrl$monitor,
               package = "quadrupen")
  ## END OF CALL
  if (ctrl$timer) {
    internal.timer <- (proc.time() - cpp.start)[3]
    external.timer <- (proc.time() - r.start)[3]
  } else {
    internal.timer <- NULL
    external.timer <- NULL
  }

  ## ======================================================                                                                                                                                                  ## BUILDING THE QUADRUPEN OBJECT
  out$converge[out$converge == 0] <- "converged"
  out$converge[out$converge == 1] <- "max # of iterate reached"
  out$converge[out$converge == 2] <- "max # of feature reached"
  monitoring  <- list()
  monitoring$it.active <- c(out$it.active)
  monitoring$it.optim  <- c(out$it.optim)
  monitoring$max.grad  <- c(out$max.grd)
  monitoring$status    <- c(out$converge)
  monitoring$nbr.in    <- c(out$nbr.in)
  monitoring$pensteps.timer <- c(out$timing)
  monitoring$external.timer <- external.timer
  monitoring$internal.timer <- internal.timer
  monitoring$dist.to.opt    <- c(out$delta.hat)
  monitoring$dist.to.str    <- c(out$delta.star)
  dim.names <- list()
  dim.names[[1]] <- round(c(out$lambda1),3)
  dim.names[[2]] <- 1:p ## ifelse(is.null(colnames(x)),1:p,colnames(x))
  coefficients <- sparseMatrix(i = out$iA+1,
                               j = out$jA+1,
                               x = c(out$nzeros),
                               dims = c(length(out$lambda1),p),
                               dimnames = dim.names)

  ## FITTED VALUES AND RESIDUALS...
  if (intercept) {
    mu <- Matrix(rep(t(c(out$mu)),n), nrow=n, byrow=TRUE)
    fitted <- mu + x %*% t(coefficients)
  } else {
    fitted <- x %*% t(coefficients)
  }
  residuals <- sweep(matrix(rep(y, length(out$lambda1)), nrow=n) - fitted, 2L, 1/sqrt(weights), "/", check.margin = FALSE)
  
  return(new("quadrupen",
             coefficients = coefficients   ,
             intercept    = intercept      ,
             normalize    = normalize      ,
             mu           = c(out$mu)      ,
             meanx        = input$xbar     ,
             normx        = input$normx    ,
             fitted       = fitted         ,
             residuals    = residuals      ,
             penscale     = penscale       ,
             penalty      = "elastic net"  ,
             naive        = naive          ,
             lambda1      = c(out$lambda1) ,
             lambda2      = lambda2        ,
             struct       = struct         ,
             weights      = weights        ,
             monitoring   = monitoring     ,
             control      = ctrl))
  
}


standardize <- function(x,y,weights,intercept,normalize,zero=.Machine$double.eps) {

  ## ============================================
  ## INTERCEPT AND NORMALIZATION TREATMENT
  if (intercept) {
    if (any(weights != 1)) {
      xbar <- apply(x, 2, weighted.mean, weights)
      ybar <- weighted.mean(y, weights)
    } else {
      xbar <- colMeans(x)
      ybar <- mean(y)
    }
  } else {
    xbar <- rep(0,ncol(x))
    ybar <- 0
  }
  
  ## ============================================
  ## NORMALIZATION
  if (normalize) {
    if (any(weights != 1)) {
      normx <- sqrt(drop(crossprod(x^2,weights)-xbar^2))
    } else {
      normx <- sqrt(drop(colSums(x^2)- xbar^2))
    }
    if (any(normx < zero)) {
      warning("A predictor has no signal: you should remove it.")
      normx[abs(normx) < zero] <- 1 ## dirty way to handle 0/0
    }
    ## xbar is scaled to handle internaly the centering of X for sparsity purpose
    xbar <- xbar/normx
  } else {
    if (any(weights != 1)) {
      normx <- sqrt(drop(crossprod(x^2, weights)))
    } else {
      normx <- rep(1,ncol(x))
    }
  }
  normy <- sqrt(sum(y^2 * weights))
  
  ## Building the sparsely encoded design matrix
  if (inherits(x, "sparseMatrix")) {
    xs    <- as(sweep(sqrt(weights)*x, 2L, normx, "/", check.margin = FALSE), "dgCMatrix")
    x     <- list(Xi = xs@i, Xj = xs@p, Xnp = diff(xs@p), Xx = xs@x)
    xty   <- drop(crossprod(sqrt(weights)*(y-ybar),scale(xs,xbar,FALSE)))
  } else {
    x     <- list(Xx = as.matrix(sweep(sqrt(weights)*x, 2L, normx, "/", check.margin = FALSE)))
    xty   <- drop(crossprod(sqrt(weights)*(y-ybar),scale(x$Xx,xbar,FALSE)))
  }
  
  return(list(xbar=xbar, ybar=ybar, normx=normx, normy=normy, xty=xty, x=x))
}

get.lambda1 <- function(xty,penscale,nlambda1,min.ratio) {
  ## ======================================================
  ## GENERATE A GRID OF PENALTY IF NONE HAS BEEN PROVIDED
  lmax <- max(abs(xty/penscale))
  return(10^seq(log10(lmax), log10(min.ratio*lmax), len=nlambda1))
}
