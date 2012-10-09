##' Class "quadrupen"
##'
##' Class of object returned by any fitting function of the
##' \pkg{quadrupen} package, e.g.; \code{elastic.net}.
##'
##' @section Slots: \describe{
##'
##' \item{\code{coefficients}:}{Matrix (class \code{"Matrix"}) of
##' coefficients with respect to the original input (even when
##' normalized). The number of rows of \code{coefficients} corresponds
##' the length of \code{lambda1}.}
##' 
##' \item{\code{intercept}:}{logical; indicates if an intercept has
##'  been included to the model.}
##'
##' \item{\code{normalize}:}{logical; indicates if the matrix of
##'  features has been normalized for the fitting.}
##' 
##' \item{\code{mu}:}{A vector (class \code{"numeric"})
##' containing the successive value of the (unpenalized) intercept.
##' Equals to zero if \code{intercept} has been set to \code{FALSE}.}
##'
##' \item{\code{meanx}:}{Vector (class \code{"numeric"}) containing the
##' mean of the design ' matrix. Vector of zero if \code{intercept} has been
##'  set to \code{FALSE}.}
##'
##' \item{\code{normx}:}{Vector (class \code{"numeric"}) containing the
##' square root of the sum of squares of each columns of the design
##' matrix. Vector of one if normalized has been set to \code{FALSE}.}
##' 
##' \item{\code{penscale}:}{Vector \code{"numeric"} with real positive
##' values that weight the l1-penalty of each feature.}
##' 
##' \item{\code{penalty}:}{Object of class \code{"character"}
##' indicating the penalizer used (only "elastic-net" available for
##' the moment).}
##' 
##' \item{\code{naive}:}{logical; indicates whether the 'naive' or the
##' 'classic' elastic-net has been fitted.}
##'
##' \item{\code{lambda1}:}{Vector (class \code{"numeric"}) of
##' l1-penalty levels for which the model has eventually been fitted.}
##' 
##' \item{\code{lambda2}:}{Scalar (class \code{"numeric"}) for the
##' amount of ridge-like penalty.}
##' 
##' \item{\code{struct}:}{\code{"Matrix"} structuring the
##' coefficients, possibly sparsely encoded (MUST be positive
##' definite).}
##' 
##' \item{\code{control}:}{Object of class \code{"list"}.}
##' 
##' \item{\code{monitoring}:}{List (class \code{"list"}) which
##' contains various indicator dealing with the optimization process.}
##' 
##' \item{\code{residuals}:}{Vector (class \code{"numeric"}) of (weighted)
##' residuals.}
##' 
##' \item{\code{fitted}:}{Vector (class \code{"numeric"}) of fitted values.}
##' 
##' \item{\code{weights}:}{Vector of observation weights used.}
##' }
##'
##' @section Methods:
##' This class comes with the usual \code{predict(object, newx, ...)},
##' \code{fitted(object, ...)}, \code{residuals(object, ...)},
##' \code{print(object, ...)}, \code{show(object)} and
##' \code{deviance(object, ...)} generic (undocumented) methods.
##'
##' Two specific and documented methods are also available: 
##'  
##' \describe{
##'    \item{plot}{plots a \code{quadrupen} object}
##'    \item{crossval}{cross-validates a \code{quadrupen} object}
##'  }
##' 
##' @aliases fitted,quadrupen-method predict,quadrupen-method
##' deviance,quadrupen-method print,quadrupen-method
##' show,quadrupen-method residuals,quadrupen-method
##'
##' @docType class
##' 
##' @keywords class
##'
##' @seealso See also \code{\link{plot.quadrupen}}, and
##' \code{\link{crossval}}.
##' 
##' @name quadrupen-class
##' @rdname quadrupen-class
##'
##' @exportClass quadrupen
##' @exportMethod fitted
##' @exportMethod residuals
##' @exportMethod predict
##' @exportMethod deviance
##' @exportMethod print
##' @exportMethod show
##'
##' @importFrom stats fitted predict residuals deviance
##' 
setClass("quadrupen",
  representation = representation(
     coefficients  = "Matrix"   ,
     intercept     = "logical"  ,
     normalize     = "logical"  ,
     mu            = "numeric"  ,
     meanx         = "numeric"  ,
     normx         = "numeric"  ,
     weights       = "numeric"  ,
     fitted        = "Matrix"   ,
     residuals     = "Matrix"   ,
     penscale      = "numeric"  ,
     penalty       = "character",
     naive         = "logical"  ,
     lambda1       = "numeric"  ,
     lambda2       = "numeric"  ,
     struct        = "Matrix"   ,
     control       = "list"     ,
     monitoring    = "list")
)

setMethod("print", "quadrupen", definition =
   function(x, ...) {
     ncoef <- ncol(x@coefficients)
     if (is.null(ncoef)) {ncoef <- ncol(x@coefficients)}
     if (x@naive) {
       cat("Linear regression with", x@penalty, "penalizer, no rescaling of the coefficients (naive).\n")
     } else {
       cat("Linear regression with", x@penalty, "penalizer, coefficients rescaled by (1+lambda2).\n")
     }
     if (any(x@intercept != 0)) {
       cat("- number of coefficients:", ncoef,"+ intercept\n")
     } else {
       cat("- number of coefficients:", ncoef,"(no intercept)\n")
     }
     cat("- penalty parameter lambda1:", length(x@lambda1), "points from",
         format(max(x@lambda1), digits = 3),"to",
         format(min(x@lambda1), digits = 3),"\n")
     cat("- penalty parameter lambda2:", x@lambda2)
     cat("\n")
     invisible(x)
   }
)

setMethod("show", "quadrupen", definition = 
   function(object) {print(object)}
)

##' Plot method for a quadrupen object
##'
##' Produce a plot of the solution path of a \code{quadrupen} fit.
##'
##' @usage plot.quadrupen(x, y, xvar = "lambda",
##'         main = paste(slot(x, "penalty")," path", sep=""),
##'         log.scale = TRUE, standardize=TRUE, reverse=FALSE,
##'         labels = NULL, plot = TRUE, ...)
##' @param x output of a fitting procedure of the \pkg{quadrupen}
##' package (only \code{\link{elastic.net}} for the moment). Must be
##' of class \code{quadrupen}.
##' @param y used for S4 compatibility.
##' @param xvar variable to plot on the X-axis: either \code{"lambda"}
##' (l1-penalty level) or \code{"fraction"} (l1-norm of thje
##' coefficients). Default is set to \code{"lambda"}.
##' @param main the main title. Default is set to the model name followed
##' by what is on the Y-axis.
##' @param log.scale logical; indicates if a log-scale should be used
##' when \code{xvar="lambda"}.
##' @param standardize logical; standardize the coefficient before
##' plotting. Default is \code{TRUE}.
##' @param reverse logical; should the X-axis be reversed when
##' \code{xvar="lambda"}? Default is \code{FALSE}.
##' @param label vector indicating the names associated to the plotted
##' variables. When specified, a legend is drawn in order to identify
##' each variable. Only relevant when the number of predictor is
##' small. Remind that the intercept does not count. Default is
##' \code{NULL}.
##' @param plot logical; indicates if the graph should be plotted on
##' call. Default is \code{TRUE}.
##'
##' @return a ggplot2 object which can be plotted via the \code{print}
##' method.
##' @seealso \code{\linkS4class{quadrupen}}.
##' 
##' @name plot,quadrupen-method
##' @aliases plot,quadrupen-method
##' @aliases plot.quadrupen
##' @docType methods
##' @rdname plot.quadrupen
##'
##' @examples \dontrun{
##' rm(list=ls())
##' library(quadrupen)
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
##' Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
##' diag(Sigma) <- 1
##' ## This gives a great advantage to the elastic-net
##' ## for support recovery; have a look at image(Sigma)
##' n <- 100
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##' lasso <- elastic.net(x,y, lambda2=0)
##' enet  <- elastic.net(x,y, lambda2=10)
##' 
##' ## Plot the Lasso path
##' plot(lasso, main="Lasso solution path")
##' ## Plot the Elastic-net path
##' plot(enet, main = "Elastic-net solution path")
##' ## Plot the Elastic-net path (fraction on X-axis, unstandardized coefficient)
##' plot(enet, standardize=FALSE, xvar="fraction")
##' }
##' 
##' @importFrom graphics plot
##' @exportMethod plot
##' @import ggplot2 scales grid
##' @export
setMethod("plot", "quadrupen", definition =
   function(x, y, xvar = "lambda",
            main = paste(slot(x, "penalty")," path", sep=""),
            log.scale = TRUE, standardize=TRUE, reverse=FALSE,
            labels = NULL, plot = TRUE, ...) {
     
     if (length(x@lambda1) == 1) {
       stop("Not available when length(lambda1) == 1")
     }
     
     nzeros <- which(colSums(x@coefficients) != 0) 
     if (length(nzeros) == 0) {
       stop("Nothing to plot: all coefficients are zero.")
     }

     beta  <- as.matrix(x@coefficients[, nzeros])
     rownames(beta) <- NULL ## avoid warning message in ggplot2
     if (standardize) {
       beta  <- scale(beta, FALSE, 1/x@normx[nzeros])
     }
     
     xv <- switch(xvar,"fraction" = apply(abs(beta),1,sum)/max(apply(abs(beta),1,sum)), x@lambda1)
     if (log.scale & xvar=="lambda") {
       xv <- log10(xv)
     }          
     
     data.coef <- melt(data.frame(xvar=xv, beta=beta),id="xvar")
     if (is.null(labels)) {
       data.coef$labels <- factor(rep(nzeros, each=length(xv)))
     } else {
       if (sum(is.na(labels[nzeros]))>0 ) {
         labels <- NULL
         warning("The number of label is wrong, ignoring them.") 
         data.coef$labels <- factor(rep(nzeros, each=length(xv))) 
       } else {
         data.coef$labels <- factor(rep(labels[nzeros], each=length(xv)))
       }
     }
     colnames(data.coef) <- c("xvar","var","coef", "variables")
     
     d <- ggplot(data.coef,aes(x=xvar,y=coefficients, colour=variables, group=var)) +
       geom_line(aes(x=xvar,y=coef)) +
         labs(x=ifelse(xvar=="fraction",expression(paste("|",beta[lambda[1]],"|",{}[1]/max[lambda[1]],"|",beta[lambda[1]],"|",{}[1],sep="")),
                ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1]))),
              y=ifelse(standardize, "standardized coefficients","coefficients")) + ggtitle(main) +
           geom_hline(yintercept=0, alpha=0.5, linetype="dotted")
     if (xvar=="lambda" & reverse==TRUE) {
       d <- d + scale_x_reverse()
     }
     if (is.null(labels)) {
       d <- d + theme(legend.position="none")
     } else {
       if (length(labels[nzeros]) != length(nzeros)) {
         d <- d + theme(legend.position="none")         
       }
     }
     if (plot) {print(d)}
     return(d)
     
   }
)

setMethod("fitted", "quadrupen", definition =
   function(object, ...) {
     return(object@fitted)
   }
)

setMethod("predict", "quadrupen", definition =
   function (object, newx=NULL, ...)  {
     if (is.null(newx)) {
       return(object@fitted)
     } else {
       if (length(object@intercept) != 1) {
         intercept <- matrix(rep(t(object@intercept), nrow(newx)),
                             nrow=nrow(newx), byrow=TRUE)
         return(intercept + newx %*% t(object@coefficients))
       } else {
         return(newx %*% t(object@coefficients))
       }
     }
   }
)

setMethod("residuals", "quadrupen", definition =
   function(object, newx=NULL, newy=NULL, ...) {
     if (is.null(newx) | is.null(newy)) {
       return(object@residuals)
     } else {
       n <- length(object@lambda1)
       return(matrix(rep(newy, n), ncol=n) - predict(object, newx))
     }
   }
)

setMethod("deviance", "quadrupen", definition =
   function(object, newx=NULL, newy=NULL, ...) {
     dev <- apply(residuals(object, newx, newy)^2,2,sum)
     return(dev)     
   }
)
