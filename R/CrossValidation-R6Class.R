#' Class CrossValidation
#'
#' Class of object returned by the [`QuadrupenFit$cross_validate()`][QuadrupenFit] method or the
#' [`cross_validate()`] function.  Owns [print()] and [plot()] methods.
#'
#' @importFrom dplyr filter
#'
#' @export
CrossValidation <- R6::R6Class(
  classname = "CrossValidation",
  private = list(
    cv_min = NA,
    cv_1se = NA,
    folds_ = NA,
    value  = NA
  ),
  active = list(
    #' @field data a data frame containing the mean
    #' cross-validated error and its associated standard error for each
    #' values of \code{lambda1} and \code{lambda2}.
    data        = function(value) {private$value},
    #' @field folds list of \code{K} vectors indicating the folds
    #' used for cross-validation.
    folds       = function(value) {private$folds_},
    #' @field lambda1 vector of \eqn{\lambda_1}{lambda1}
    #' (\eqn{\ell_1}{l1} or \eqn{\ell_\infty}{infinity} penalty levels)
    #' for which each cross-validation has been performed.
    lambda1     = function(value) unique(self$data$lambda1),
    #' @field lambda2 vector (or scalar) of \eqn{\ell_2}{l2}-penalty levels for
    #' which each cross-validation has been performed.
    lambda2     = function(value) unique(self$data$lambda2),
    #' @field lambda1_min level of \eqn{\lambda_1}{lambda1} that minimizes the
    #' error estimated by cross-validation.
    lambda1_min = function(value) max(private$cv_min$lambda1),
    #' @field lambda2_min level of \eqn{\lambda_2}{lambda2} that minimizes the
    #' error estimated by cross-validation.
    lambda2_min = function(value) max(private$cv_min$lambda2),
    #' @field lambda1_1se largest level of \eqn{\lambda_1}{lambda1} such as
    #' the cross-validated error is within 1 standard error of the
    #' minimum.
    lambda1_1se = function(value) max(private$cv_1se$lambda1),
    #' @field lambda2_1se largest level of \eqn{\lambda_2}{lambda2}
    #' such that the cross-validated error is within 1 standard error of the
    #' minimum (only relevant for ridge regression).
    lambda2_1se = function(value) max(private$cv_1se$lambda2)
  ),
  public = list(
    #' @description Constructor for a [CrossValidation] object
    #' Should be called internally by an object [`QuadrupenFit$cross_validate()`][QuadrupenFit]
    #' @param cv_error data frame storing output of a cv job
    #' @param folds  list of K folds used for cross-validation
    initialize =
      function(cv_error, folds) {
        private$value  <- cv_error
        private$folds_ <- folds
        private$cv_min <- cv_error |> dplyr::filter(mean <= min(mean))
        private$cv_1se <- cv_error |>
          dplyr::filter(lambda2 == self$lambda2_min) |>
          dplyr::filter(mean < min(mean + serr) + 1e-5)
    },
    #' @description User friendly print method
    show = function() {
      cat(length(self$folds), "fold cross-validation job for a Quadrupen fit.\n")
      cat("- main penalty parameter:", length(self$lambda1), "points from",
          format(max(self$lambda1), digits = 3),"to",
          format(min(self$lambda1), digits = 3),"\n")
      if (length(self$lambda2) > 1) {
        cat("- structuring penalty parameter:", length(self$lambda2), "points from",
            format(max(self$lambda2), digits = 3),"to",
            format(min(self$lambda2), digits = 3),"\n")
      } else {
        cat("- structuring penalty parameter:", self$lambda2, "\n")
      }
      invisible(self)
    },
    #' @description User friendly print method
    print = function() { self$show() },
    #' @description Plot 1-dimensional cross-validation
    #' @param log_scale logical, should a log-scale be used for the x-axis
    #' @param se logical, should confidence band be displayed (TRUE by default)
    #' @param title graph title
    #' @return a ggplot object
    plotCV_1D = function(log_scale=TRUE, title = "Cross-validation error", se = TRUE) {
      ## SIMPLE CROSS-VALIDATION GRAPH
      if (length(self$lambda1) > length(self$lambda2)) {
        d <- ggplot(self$data, aes(x=.data$lambda1,y=.data$mean)) + theme_bw()
        xlab <- ifelse(log_scale,expression(log[10](lambda[1])),expression(lambda[1]))
        lambda <- data.frame(xval=c(self$lambda1_min,self$lambda1_1se),
                             lambda.choice=factor(c("min. MSE","1-se rule")))

      } else {
        ## ridge or not (meaning working on lambda1 or lambda2)
        d <- ggplot(self$data, aes(x=.data$lambda2,y=.data$mean))
        xlab <- ifelse(log_scale,expression(log[10](lambda[2])),expression(lambda[2]))
        lambda <- data.frame(xval=c(self$lambda2_min,self$lambda2_1se),
                             lambda.choice=factor(c("min. MSE","1-se rule")))
      }
      d <- d + xlab(xlab) + ylab("Mean square error") +
        geom_point(alpha=0.3) +
        geom_smooth(aes(ymin=.data$mean-.data$serr, ymax=.data$mean+.data$serr, group=as.factor(.data$lambda2),
                        colour=as.factor(.data$lambda2)), se=se, data=self$data, alpha=0.2, stat="identity")
      if (log_scale) {
        d <- d + scale_x_log10() + annotation_logticks(sides="b")
      }
      d <- d + ggtitle(title) +
        geom_vline(data=lambda, aes(xintercept=.data$xval, linetype=.data$lambda.choice),
                   alpha=0.2, show.legend = TRUE) +
          scale_color_discrete(labels = round(self$lambda2, 3)) +
        labs(colour = expression(explored~lambda[2]), linetype = expression(lambda[1]~choice))

      d
    },
    #' @description Plot 2-dimensional cross-validation output (grid lambda1 x lambda2)
    #' @param title graph title
    #' @return a \pkg{ggplot2} object
    plotCV_2D = function(title = "Cross-validation error") {
      d <- ggplot(data=self$data, aes(x=.data$lambda1, y=.data$lambda2, z=.data$mean))
      d <- d + geom_tile(aes(fill=.data$mean)) + stat_contour(linewidth=0.2, binwidth=diff(range(self$data$mean))/10) + ggtitle(title)
      d <- d + scale_x_continuous(trans=log10_trans()) + xlab(expression(log[10](lambda[1])))
      d <- d + scale_y_continuous(trans=log10_trans()) + ylab(expression(log[10](lambda[2])))
      d <- d + annotation_logticks() + theme_bw()
      i_1se <- which(self$data$mean - self$data$serr <= min(self$data$mean))
      d <- d + stat_contour(alpha=0.5, colour="#CCCCCC", linewidth=0.65, breaks=quantile(self$data$mean[i_1se], probs=seq(0,1,len=6)))
      d
    },
    #' @description Plot cross-validation job by choosing the most appropriate output (1D- or 2D)
    #' @param log_scale logical, should a log-scale be used for the x-axis
    #' @param title graph title
    #' @return a ggplot object
    plot = function(log_scale=TRUE, title = "Cross-validation error") {
      if (length(self$lambda1) > 1 & length(self$lambda2) > 1) {
        d <- self$plotCV_2D(title)
      } else {
        d <- self$plotCV_1D(log_scale, title)
      }
      d
    }
  )
)
