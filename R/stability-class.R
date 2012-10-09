##' Class "stability.path"
##'
##' Class of object returned by the \code{stability} function, with
##' methods \code{print}, \code{show} and \code{plot}.
##'
##' @section Slots: \describe{
##'
##' \item{\code{probabilities}:}{a \code{Matrix} object containing the
##' estimated selection of probabilities along the l1-penalty
##' levels.}
##' \item{\code{penalty}:}{Object of class \code{"character"}
##' indicating the penalizer used.}
##' \item{\code{naive}:}{logical indicating whether rescaling of the
##' coefficients has been performed regarding the l2-penalty.}
##' \item{\code{lambda1}:}{a vector with the l1-penalty levels.}
##' \item{\code{lambda2}:}{a scalar with the l2-penalty level.}
##' }
##' @aliases print,stability.path-method show,stability.path-method
##'
##' @docType class
##' 
##' @keywords class
##'
##' @seealso See also \code{\link{plot.stability.path}}, and
##' \code{\link{stability}}.
##' 
##' @name stability.path-class
##' 
##' @rdname stability.path-class
##'
##' @exportClass stability.path
##' @exportMethod print
##' @exportMethod show
##'
setClass("stability.path",
  representation = representation(
     probabilities = "Matrix"   ,
     penalty       = "character",
     naive         = "logical"  ,
     lambda1       = "numeric"  ,
     lambda2       = "numeric")
)

setMethod("print", "stability.path", definition =
   function(x, ...) {
     if (x@naive) {
       cat("Stability path for", x@penalty, "penalizer, no rescaling of the coefficients (naive).\n")
     } else {
       cat("Stability path for", x@penalty, "penalizer, coefficients rescaled by (1+lambda2).\n")
     }
     cat("- penalty parameter lambda1:", length(x@lambda1), "points from",
         format(max(x@lambda1), digits = 3),"to",
         format(min(x@lambda1), digits = 3),"\n")
     cat("- penalty parameter lambda2:", x@lambda2)
     cat("\n")
     invisible(x)
   }
)

setMethod("show", "stability.path", definition = 
   function(object) {print(object)}
)

##' Plot method for \code{stability.path}.
##'
##' Produce a plot of the stability path obtained by stability
##' selection. Display the average number of selected variables to
##' controlled the per-family error rate at a given cutoff
##' probability.
##'
##' @usage plot.stability.path(x, y, xvar = "ave.sel", annot=TRUE,
##'          main = paste("Stability path of an ", slot(x, "penalty")," regularizer", sep=""),
##'          log.scale = TRUE,  labels = rep("unknown status",p), plot = TRUE, cutoff=0.75, PFER=2, ...)
##' @param x output of a \code{crossval} run (must be of class
##' \code{stability.path}).
##' @param y used for S4 compatibility.
##' @param main main title.
##' @param xvar variable to plot on the X-axis: either
##' \code{"ave.sel"} (average number of selected features),
##' \code{"lambda"} (l1-penalty level) or \code{"fraction"} (fraction
##' of l1-penalty level applied). Default is set to \code{"ave.sel"}.
##' @param log.scale logical; indicates if a log-scale should be used
##' when \code{xvar="lambda"}. 
##' @param plot logical; indicates if the graph should be
##' plotted. Default is \code{TRUE}.
##' @param annot logical; should annotation be made on the graph
##' regarding controlled PFER. Default is \code{TRUE}.
##' @param cutoff value of the cutoff probability. See details.
##' @param PFER value of the per-family error rate to control. See
##' details.
##' @param labels an optional vector of labels for each variable in
##' the path (e.g., 'relevant'/'irrelevant').
##' @param ... used for S4 compatibility.
##' @return a list with a ggplot2 object which can be plotted via the
##' \code{print} method, and a vector of selected variables corresponding to the
##' given cutoff and PFER argument.
##' @name plot,stability.path-method
##' @aliases plot,stability.path-method
##' @aliases plot.stability.path
##' @docType methods
##' @rdname plot.stability.path
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
##' n <- 100
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##' 
##' ## Build a vector of label for true nonzeros
##' labels <- rep("irrelevant", length(beta))
##' labels[beta != 0] <- c("relevant")
##' labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))
##' 
##' ## Call to stability selection function, 200 subsampling
##' stabout <- stability(x,y, subsamples=200, lambda2=1, min.ratio=1e-2)
##' ## Build the plot an recover the selected variable for a given cutoff
##' ## and per-family error rate
##' stabpath <- plot(stab, labels=labels, cutoff=0.75, PFER=1)
##' stabpath <- plot(stab, xvar="lambda", labels=labels, cutoff=0.75, PFER=2)
##' stabpath <- plot(stab, xvar="fraction", annot=FALSE, labels=labels, cutoff=0.75, PFER=1)
##' }
##' @importFrom graphics plot
##' @exportMethod plot
##' @import ggplot2 scales grid
##' @export
setMethod("plot", "stability.path", definition =
   function(x, y, xvar = "ave.sel", annot=TRUE,
            main = paste("Stability path of an ", slot(x, "penalty")," regularizer", sep=""),
            log.scale = TRUE,  labels = rep("unknown status",p), plot = TRUE, cutoff=0.75, PFER=2, ...) {

     p <- ncol(x@probabilities)
     if (length(x@lambda1) == 1)
       stop("Not available when length(lambda1) == 1")

     if(PFER <=0)
       stop("PFER should be at least equal to 1.")

     if(cutoff <0.5 | cutoff > 1)
       stop("The cutoff is supposed to be a probability in [.5,1] ...")

     
     nzeros <- which(colSums(x@probabilities) != 0) 
     if (length(nzeros) == 0)
       stop("Nothing to plot: all probabilities are zero along the path.")

     prob  <- as.matrix(x@probabilities[, nzeros])
     rownames(prob) <- NULL

     ## estimate the average number of selected variables on the
     ## current path and pick the one controling the PFER at the
     ## desired level
     q <- cumsum(rowSums(prob))/1:nrow(prob)
     iq <- which.min(sqrt(ncol(prob) * (2*cutoff-1) * PFER) >= q)

     ## the x-axis variable
     xv <- switch(xvar,
                  "fraction" = x@lambda1/max(x@lambda1),
                  "ave.sel" = q, x@lambda1)
     if (log.scale & xvar=="lambda")
       xv <- log10(xv)

     ## Build the data frame for ggploting
     selection <- rep("unselected",ncol(prob))
     selection[prob[iq, ] > cutoff] <- "selected"
     data.coef <- melt(data.frame(xvar=xv, prob=prob),id="xvar")
     data.coef$selection <- factor(rep(selection, each=length(xv)))
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
     colnames(data.coef) <- c("xvar","var","prob","selection","variables")

     ## Build the ggplot object
     d <- ggplot(data.coef,aes(x=xvar,y=prob, linetype=variables, colour=selection, group=var)) +
       geom_line(aes(x=xvar,y=prob)) +
         labs(x=switch(xvar,
                "fraction" = expression(lambda[1]/max[lambda[1]]),
                "ave.sel" = "average number of selected variables",
                ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1]))),
              y="selection probabilities") + ggtitle(main) 
     if (xvar!="ave.sel") {
       d <- d + scale_x_reverse()
     } 
     if (is.null(labels)) {
       d <- d + theme(legend.position="none")
     } else {
       if (length(labels[nzeros]) != length(nzeros)) {
         d <- d + theme(legend.position="none")         
       }
     }

     ## Manage the annotation for the selected variables
     if (annot) {
       d <- d + annotate("rect", xmin=xv[1], xmax=xv[iq], ymin=cutoff, ymax=1, alpha=0.15)
       d <- d + annotate("text", x=c(xv[length(xv)],xv[1],xv[iq]), y=c(cutoff,1,0), hjust=c(0,0,0.25), vjust=c(0,1.1,1.1),
                         label=c(paste("pi[thr]"),paste("PFER <=",PFER),paste("hat(q)==",round(q[iq],2))),
                         parse=TRUE, size=4, alpha=0.85)
       d <- d + geom_hline(yintercept=cutoff, linetype="dashed", alpha=0.35, size=.5)
       d <- d + geom_vline(xintercept=xv[iq], linetype="dashed", alpha=0.35, size=.5)
     }
     
     if (plot) {print(d)}
     invisible(list(ggplot.object=d, selected=which(selection == "selected")))
     
   }
)

