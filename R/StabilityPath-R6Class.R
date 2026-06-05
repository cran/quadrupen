#' Class StabilityPath
#'
#' Class of object returned by the [`QuadrupenFit$cross_validate()`][QuadrupenFit] method or the
#' [`cross_validate()`] function. Owns [print()] and [plot()] methods.
#'
#' @param sel_mode a character string, either `"rank"` or
#' `"PFER"`. In the first case, the selection is based on the
#' rank of total probabilities by variables along the path: the first
#' `nvarsel` variables are selected (see below). In the second
#' case, the PFER control is used as described in Meinshausen and
#' Buhlmannn's paper. Default is `"rank"`.
#' @param nvarsel number of variables selected (only relevant when
#' `sel_mode` equals `"rank"`. Default is `floor(n/log(p))`.
#' @param cutoff value of the cutoff probability (only relevant when
#' `sel_mode` equals `"PFER"`).
#' @param PFER value of the per-family error rate to control (only
#' relevant when `sel_mode` equals `"PFER"`).
#' 
#' @export
StabilityPath <- R6::R6Class(
  classname = "StabilityPath",
  public = list(
    #' @field probabilities a `Matrix` object containing the
    #' estimated probabilities of selection along the path of solutions.
    probabilities = NA, 
    #' @field regParam a list with the levels of the regularizing parameters used
    regParam      = NA, 
    #' @field subsamples a list that contains the folds used for each subsample.
    subsamples    = NA, 
    ## model-related fields
    #' @description Constructor for a [StabilityPath] object
    #' Should be called internally by an object [`QuadrupenFit$stability()`][QuadrupenFit]
    #' @param probabilities a `Matrix` object containing the
    #' estimated probabilities of selection along the path of solutions.
    #' @param regParam a list with the levels of the regularizing parameters used
    #' @param subsamples a list that contains the folds used for each subsample.
    initialize = 
      function(probabilities, regParam, subsamples) {
        self$probabilities <- probabilities
        self$regParam      <- regParam
        self$subsamples    <- subsamples
      },
    #' @description User friendly print method
    show = function() {
      cat("Stability path with", length(self$subsamples), "resamplings.\n")
      cat("- penalty parameter lambda1:", length(self$regParam[[1]]), "points from",
          format(max(self$regParam[[1]]), digits = 3),"to",
          format(min(self$regParam[[1]]), digits = 3),"\n")
      if (!is.null(self$regParam[[2]]))
        cat("- penalty parameter lambda2:", self$regParam[[2]], "\n")
      invisible(self)
    },
    #' @description User friendly print method
    print = function() { self$show() },
    #' @description Perform variable selection based on the stability path
    selection = function(sel_mode = c("rank", "PFER"), cutoff=0.75, 
                         PFER=2, nvarsel=floor(self$nobs/log(self$nvar))) {

      n <- self$nobs
      p <- self$nvar
      
      stopifnot("PFER should be at least equal to 1." = (PFER > 0))
      stopifnot("The cutoff is supposed to be a probability in [.5,1] ..." = (cutoff >= 0.5 &  cutoff <= 1))
      stopifnot("The rank is supposed to be less than p ..." = (nvarsel <= p))
      
      sel_mode  <- match.arg(sel_mode)
      prob      <- self$nonzeroprob
      
      if (sel_mode == "PFER") {
        ## estimate the average number of selected variables on the current path
        ## and pick the one controlling the PFER at the desired level
        q    <- colSums(prob >= cutoff)
        qLim <- sqrt(PFER * (2 * cutoff - 1) * p)
        iq <- q <= qLim
        iq <- ifelse(which.min(iq) != 1,which.min(iq) - 1,ifelse(sum(iq) == 0,1,length(iq)))
        selected <- self$nonzero[which(prob[, iq] > cutoff)]
        attr(selected, "iq") <- iq
        attr(selected, "q") <- q
      } else {
        id <- order(rowSums(prob),decreasing = TRUE)[1:nvarsel]
        selected <- self$nonzero[id]
      }
      selected
    },
    #' @description Produce a plot of the stability path obtained by stability
    #' selection.
    #'
    #' @param title title title. If none given, a somewhat appropriate
    #' title is automatically generated.
    #' @param xvar variable to plot on the X-axis: either `"lambda"`
    #' (first penalty level) or `"fraction"` (fraction of the
    #' penalty level applied tune by \eqn{\lambda_1}{lambda1}). Default
    #' is `"lambda"`.
    #' @param plot logical; indicates if the graph should be
    #' plotted. Default is `TRUE`. If `FALSE`, only the
    #' \pkg{ggplot2} object is sent back.
    #' @param labels an optional vector of labels for each variable in
    #' the path (e.g., 'relevant'/'irrelevant'). See examples.
    #' @return a list with a \pkg{ggplot2} object which can be plotted
    #' via the \code{print} method, and a vector of selected variables
    #' corresponding to method of choice (`"rank"` or
    #' `"PFER"`).
    #'
    #' @examples \dontrun{
    #' ## Simulating multivariate Gaussian with blockwise correlation
    #' ## and piecewise constant vector of parameters
    #' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
    #' Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
    #' Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
    #' Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
    #' diag(Sigma) <- 1
    #' n <- 100
    #' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
    #' y <- 10 + x %*% beta + rnorm(n,0,10)
    #'
    #' ## Build a vector of label for true nonzeros
    #' labels <- rep("irrelevant", length(beta))
    #' labels[beta != 0] <- c("relevant")
    #' labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))
    #'
    #' enet <- elastic_net(x, y, lambda2 = 10, struct = solve(Sigma), minratio = 1e-2)
    #' stab <- stability(enet, n_subsamples = 200)
    #'
    #' ## Build the plot an recover the selected variable
    #' plot(stab, labels=labels)
    #' 
    #' cat("\nFalse positives for the randomized Elastic-net with stability selection: ",
    #'      sum(labels[stab$selection()] != "relevant"))
    #' cat("\nDONE.\n")
    #' }
    #' @importFrom graphics plot
    #' @import ggplot2 scales grid
    #' 
    plot = function(
      xvar = "lambda", title = "Stability path", labels = rep("unknown status",self$nvar), 
      sel_mode = c("rank", "PFER"), cutoff=0.75, PFER=2, nvarsel=min(self$nvar, floor(self$nobs/log(self$nvar)))) {

      sel_mode <- match.arg(sel_mode)
      selection <- rep("unselected",length(self$nonzero))
      selected <- self$selection(sel_mode, cutoff, PFER, nvarsel)
      selection[selected] <- "selected"
      
      stopifnot("Not available when length(lambda1) == 1" = (length(self$regParam[[1]]) > 1))
      stopifnot("Nothing to plot: all probabilities are zero along the path." = (length(self$nonzero) > 0))

      prob <- t(self$nonzeroprob)
      
      ## the x-axis variable
      xv <- switch(xvar,"fraction" = self$regParam[[1]]/max(self$regParam[[1]]), self$regParam[[1]])
      if (xvar == "lambda") xv <- log10(xv)

      ## Build the data frame for ggploting
      dplot <- data.frame(xvar=xv, prob=prob) |> 
        tidyr::pivot_longer(cols = -xvar, names_to = "var", values_to = "prob") |>
        mutate(selection = factor(rep(selection, length(xv))))
      
      # dplot$selection <- factor(rep(selection, length(xv)))
      if (is.null(labels)) {
        dplot$labels <- factor(rep(1:length(self$nonzero), length(xv)))
      } else {
        if (sum(is.na(labels[self$nonzero])) > 0 ) {
          labels <- NULL
          warning("The number of label is wrong, ignoring them.")
          dplot$labels <- factor(rep(self$nonzero, length(xv)))
        } else {
          dplot$labels <- factor(rep(labels[self$nonzero], length(xv)))
        }
      }

      ## Build the ggplot object
      d <- ggplot(dplot) + 
        aes(x=xvar,y=prob, linetype=labels, color=selection, group=var) + 
        geom_line() +
        labs(x=switch(xvar,
                      "fraction" = expression(lambda/max[lambda]),
                      expression(log[10](lambda))),
             y="selection probabilities") + ggtitle(title) + theme_bw()
      d <- d + scale_x_reverse()
      if (is.null(labels)) {
        d <- d + theme(legend.position="none")
      } else {
        if (length(labels[self$nonzero]) != length(self$nonzero)) {
          d <- d + theme(legend.position="none")
        }
      }
      
      ## Manage the annotation for the selected variables (PFER mode)
      if (sel_mode == "PFER") {
        iq <- attr(selected, "iq")
        q  <- attr(selected, "q")
        d <- d + annotate("rect", xmin=xv[1], xmax=xv[iq], ymin=cutoff, ymax=1, alpha=0.15)
        d <- d + annotate("text", x=c(xv[length(xv)],xv[1],xv[iq]), y=c(cutoff,1,0), hjust=c(0,0,0.25), vjust=c(0,1.1,1.1),
                          label=c(paste("pi[thr]"),paste("PFER <=",PFER),paste("hat(q)==",round(q[iq],2))),
                          parse=TRUE, size=4, alpha=0.85)
        d <- d + geom_hline(yintercept=cutoff, linetype="dashed", alpha=0.35, linewidth=.5)
        d <- d + geom_vline(xintercept=xv[iq], linetype="dashed", alpha=0.35, linewidth=.5)
      }
      d
    }
  ),
  active =list(
    #' @field nvar number of variables (without intercept)
    nvar = function(value) {nrow(self$probabilities)},
    #' @field nobs number of observation/sample size
    nobs = function(value) {max(unlist(self$subsamples))},
    #' @field nonzero variables with a non-null probability of selection along the stability path
    nonzero = function(value) {which(rowSums(self$probabilities) != 0)},
    #' @field nonzeroprob subset of the probabilities stability path on the nonzero variables
    nonzeroprob = function(value) {
      res  <- as.matrix(self$probabilities[self$nonzero, ])
      colnames(res) <- NULL
      res      
    }
  )
)
