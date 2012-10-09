##' Cross-validation function for quadrupen fitting methods.
##'
##' Function that computes K-fold (double) cross-validated error of a
##' \code{quadrupen} fit. If no lambda2 is provided, simple cross
##' validation on the lambda1 parameter is performed. If a vector
##' lambda2 is provided, double cross-validation is performed.
##'
##' @param method a string for the fitting procedure used for
##' cross-validation. Only \code{"elastic.net"} is available for the
##' moment
##' 
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept.
##' 
##' @param y response vector.
##'
##' @param K integer indicating the number of folds. Default is 10.
##'
##' @param folds list of \code{K} vectors that describes the folds to
##' use for the cross-validation. By default, the folds are randomly
##' sampled with the specified K. The same folds are used for each
##' values of \code{lamdba2}.
##'
##' @param lambda2 tunes the l2-penalty (ridge-like) of the fit. If
##' none is provided, the default scalar value of the corresponding
##' fitting method is used and a simple CV is performed. If a vector
##' of values is provided, double cross-validation is performed (both
##' on lambda1 and lambda2, using the same folds for each lambda2).
##'
##' @param verbose logical; indicates if the progression (the current
##' fold number) should be displayed. Default is \code{TRUE}.
##'
##' @param ... additional parameters to overwrite the defaults of the
##' \code{'method'} fitting procedure. See the corresponding
##' documentation (e.g., \code{\link{elastic.net}})
##'
##' @return An object of class "cvpen" for which a \code{plot} method
##' is available.
##'
##' @seealso \code{\linkS4class{quadrupen}}, \code{\link{plot.cvpen}}
##' and \code{\linkS4class{cvpen}}.
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
##' ## Use fewer lambda1 values by overwritting the default parameters
##' ## and cross-validate over the sequences lambda1 and lambda2
##' cvout <- crossval(x,y, lambda2=10^seq(2,-2,len=50), nlambda1=50)
##' 
##' beta.1se <- slot(cvout, "beta.1se")
##' beta.min <- slot(cvout, "beta.min")
##' 
##' cat("\nFalse positives with the minimal CV choice: ", sum(sign(beta) != sign(beta.min)))
##' cat("\nFalse positives with the 1 standard-error rule: ", sum(sign(beta) != sign(beta.1se)))
##'
##' @keywords models, regression
##' @name crossval
##' @aliases crossval
##' @rdname crossval
##'
##' @export
crossval <- function(x,
                     y,
                     method  = "elastic.net",
                     K       = 10,
                     folds   = split(sample(1:nrow(x)), rep(1:K, length=nrow(x))),
                     lambda2 = 0.01,
                     verbose = TRUE, ...) {
     
  ## =============================================================
  ## INITIALIZATION & PARAMETERS RECOVERY
  user <- list(...)
  defs <- default.args(method,nrow(x),ncol(x),user)
  args <- modifyList(defs, user)
  ## Compute a grid of lambda1 (the smae for each fold)
  if (is.null(args$lambda1)) {
    input <- standardize(x,y,args$weights,args$intercept,args$normalize)
    args$lambda1 <- get.lambda1(input$xty,args$penscale,args$nlambda1,args$min.ratio)
    rm(input)
  }

  ## =============================================================
  if (length(lambda2) > 1) {
    ## DOUBLE CROSS-VALIDATION WORK
    if (verbose){
      cat("\nDOUBLE CROSS-VALIDATION\n\n")
      cat(length(folds),"-fold CV on the lambda1 grid for each lambda2\n", sep="")
    }
    cv <- sapply(1:length(lambda2), function(i) {
      if(verbose){
        cat(round(lambda2[i],3),"\t")
        if (i %% 5 == 0) {cat("\n")}
      }
      simple.cv(folds, x, y, method, args, lambda2[i])
    }, simplify=FALSE)
    if(verbose){cat("\n")}
    
    ## Recovering the best lambda1 and lambda2
    lambda1.cv <- sapply(1:length(lambda2), function(j) {       
      cv.min <- min(cv[[j]]$mean)
      lambda1.min   <- max(args$lambda1[cv[[j]]$mean <= cv.min], na.rm=TRUE)
      lambda1.1se   <- max(args$lambda1[cv[[j]]$mean <(cv[[j]]$mean+cv[[j]]$serr+1e-5)[match(lambda1.min,cv[[j]]$lambda1)]], na.rm=TRUE)
      return(c(cv.min, lambda1.min, lambda1.1se))
    })
    ind.lb2.min <- which.min(lambda1.cv[1,])
    lambda2.min <- lambda2[ind.lb2.min]
    lambda1.min <- lambda1.cv[2, ind.lb2.min]
    lambda1.1se <- lambda1.cv[3, ind.lb2.min]
    
    ## formatting cv.error for ggplot
    cv <- data.frame(do.call(rbind,cv),lambda2=rep(lambda2,sapply(cv, function(x) length(x$lambda1))))
  } else {
    ## SIMPLE CROSS-VALIDATION WORK
    if (verbose) {
      cat("\nSIMPLE CROSS-VALIDATION\n")
      cat(length(folds),"-fold CV on the lambda1 grid, lambda2 is fixed.\n", sep="")
    }
    cv <- simple.cv(folds, x, y, method, args, lambda2) 
    
    ## Recovering the best lambda1 and lambda2
    lambda1.min <- max(cv$lambda1[cv$mean <= min(cv$mean)], na.rm=TRUE)
    lambda1.1se <- max(cv$lambda1[cv$mean <(cv$mean+cv$serr+1e-5)[match(lambda1.min,cv$lambda1)]], na.rm=TRUE)
    lambda2.min <- lambda2
  }
  
  ## Apply the fitting procedure with these best lambda2 parameter
  args$lambda2 <- lambda2.min
  best.fit <- do.call(method, c(list(x=x,y=y),args))
  
  ## Finally recover the CV choice (minimum and 1-se rule)
  beta.min <- best.fit@coefficients[match(lambda1.min, args$lambda1),]
  beta.1se <- best.fit@coefficients[match(lambda1.1se, args$lambda1),]
  
  return(new("cvpen",
             lambda1     = args$lambda1,
             lambda1.min = lambda1.min,                
             lambda1.1se = lambda1.1se,
             lambda2     = lambda2,
             lambda2.min = lambda2.min,
             cv.error    = cv,
             folds       = folds,
             beta.min    = beta.min,
             beta.1se    = beta.1se))          
}

simple.cv <- function(folds, x, y, method, args, lambda2) {

  K <- length(folds)
  n <- length(y)
  y.hat <- matrix(NA, n, length(args$lambda1))

  ## overwrite irrelevant arguments
  args$control$verbose  <- FALSE
  args$lambda2  <- lambda2
  args$max.feat <- ncol(x)
  ## keep trace of the weights
  weights <- args$weights
  
  ## Multicore approach
  one.fold <- function(k) {
    omit <- folds[[k]]
    args$weights <- weights[-omit]
    fit <- do.call(method, c(list(x=x[-omit, ], y=y[-omit]), args))
    fold.err <- sweep(matrix(predict(fit,matrix(x[omit,], nrow=length(omit))), nrow=length(omit)), 1L, y[omit], check.margin = FALSE)^2
  }
  ## turn a list to matrix
  err  <- do.call(rbind,mclapply(1:K, one.fold))
  ## efficient computation of means and the standard error
  mean <- colMeans(err)
  serr <- sqrt((colSums(sweep(err, 2L, mean, check.margin = FALSE)^2)/(n-1)) /K)

  return(data.frame(mean=mean, serr=serr, lambda1=args$lambda1))
}

default.args <- function(method,n,p,user) {
  lambda2 <- ifelse(is.null(user$lambda2),0.01,user$lambda2)
  return(switch(method, "elastic.net" = list(
                          lambda1   = NULL,
                          lambda2   = 0.01,
                          penscale  = rep(1,p),
                          struct    = 1*bandSparse(p, k=0),
                          weights   = rep(1,n),
                          intercept = TRUE,
                          normalize = TRUE,
                          naive     = FALSE,
                          nlambda1  = ifelse(is.null(user$lambda1),100,length(user$lambda1)),
                          min.ratio = ifelse(n<p,0.01,5e-3),
                          max.feat  = as.integer(ifelse(lambda2<1e-2,min(n,p),min(3*n,p))),
                          control   = list(),
                          checkargs = FALSE)))
}
