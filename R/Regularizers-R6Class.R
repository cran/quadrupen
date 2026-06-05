#' Class "SparseFit"
#'
#' Class of object returned by the fitting function [elastic_net()]. Inherits fields
#' and methods of [QuadrupenFit]
#'
#' @seealso [QuadrupenFit], [elastic_net()]
#'
#' @export
#'
SparseFit <- R6::R6Class(
  classname = "SparseFit",
  inherit = QuadrupenFit,
  private  = list(type_ = NA),
  active  = list(
    #' @field lambda1 vector of tuning parameters for the l1 penalty
    lambda1 = function(value) private$tuning[[1]],
    #' @field lambda2 vector of tuning parameters for the l2 penalty
    lambda2 = function(value) private$tuning[[2]],
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) {
      ridge  <- ifelse(self$is_l2_regularized, "plus L2-regularization", "")
      sparse  <- switch(private$type_,
                       "l1" = "L1",
                       "mcp" = "MCP",
                       "scad" = "SCAD"
      )
      paste(sparse, ridge)
    },
    #' @field type string the type of group-wise regularization applied
    type = function(value) private$type_,
    #' @field unbiasing_tuning unbiasing coefficient of the MCP or SCAD penalties
    unbiasing_tuning = function(value) {
      if (missing(value))
        return(private$tuning$eta)
      else private$tuning$eta <- value
    }
  ),
  public  = list(
    #' @description Initialize a [`SparseFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, type, regParam) {
      super$initialize(data, intercept, regParam)
      private$type_  <- type
      if (data$sparse_encoding) {
        private$optimizer <-  switch(private$type_,
            "mcp" = mcp_sparse_cpp, "scad" = scad_sparse_cpp, elastic_net_sparse_cpp)
      } else {
        private$optimizer <- switch(private$type_,
            "mcp" = mcp_dense_cpp, "scad" = scad_dense_cpp, elastic_net_dense_cpp)
      }
    }
  )
)

#' Class "BoundedRegression"
#'
#' Class of object returned by the fitting function [bounded_reg()]. Inherits fields
#' and methods of [QuadrupenFit].
#'
#' @seealso [QuadrupenFit], [bounded_reg()]
#'
#' @export
BoundedRegressionFit <- R6::R6Class(
  classname = "BoundedRegressionFit",
  inherit = QuadrupenFit,
  active  = list(
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) "Bounded (LINF)",
    #' @field lambdainf vector of tuning parameters for the linf penalty
    lambdainf = function(value) private$tuning[[1]],
    #' @field lambda2 vector of tuning parameters for the l2 penalty
    lambda2 = function(value) private$tuning[[2]]
  ),
  #' @description Initialize a [`BoundedRegressionFit`] model
  #' @param data a [`DataModel`] object
  #' @param intercept a logical; should an intercept be included in the mode?
  #' @param regParam a list with two elements, a vector and a scalar, for the regularization
  public  = list(
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- bounded_regression_cpp
    }
  )
)


#' Class "RidgeRegressionFit"
#'
#' Class of object returned by the fitting function [ridge()]. Inherits fields
#' and methods of [QuadrupenFit].
#'
#' @seealso [QuadrupenFit], [ridge()]
#'
#' @export
RidgeRegressionFit <- R6::R6Class(
  classname = "RidgeRegressionFit",
  inherit = QuadrupenFit,
  active  = list(
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) "Ridge",
    #' @field lambda2 vector of tuning parameters for the l2 penalty
    lambda2 = function(value) private$tuning[[1]]
  ),
  public  = list(
    #' @description Initialize a [`RidgeRegressionFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- ridge_cpp
    }
  )
)


#' Class "FusedLassoFit"
#'
#' Class of object returned by the fitting function [fused_lasso()]. Inherits fields
#' and methods of [QuadrupenFit].
#'
#' @seealso [QuadrupenFit], [bounded_reg()]
#'
#' @export
FusedLassoFit <- R6::R6Class(
  classname = "FusedLassoFit",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(
    penalty = function(value) "Fused-LASSO",
      #' @field lambda1 vector of tuning parameters for the l1 penalty
    lambda1 = function(value) private$tuning[[1]],
    #' @field lambda2 vector of tuning parameters for the fusion penalty
    lambda2 = function(value) private$tuning[[2]]
  ),
  #' @description Initialize a [`FusedLassoFit`] model
  #' @param data a [`DataModel`] object
  #' @param intercept a logical; should an intercept be included in the mode?
  #' @param regParam a list with two elements, a vector and a scalar, for the regularization
  public  = list(
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- FusedLasso_cpp
    }
  )
)

#' Class "SparseGroupFit"
#'
#' Class of object returned by the fitting function [group_sparse_lm()]. Inherits fields
#' and methods of [QuadrupenFit]
#'
#' @seealso [QuadrupenFit], [group_sparse_lm()]
#'
#' @export
#'
SparseGroupFit <- R6::R6Class(
  classname = "SparseGroupFit",
  inherit = QuadrupenFit,
  private  = list(group_ = NA, type_ = NA),
  active  = list(
    #' @field lambda1 vector of tuning parameters for the l1 group penalty
    lambda1 = function(value) private$tuning[[1]],
    #' @field lambda2 vector of tuning parameters for the l2 penalty
    lambda2 = function(value) private$tuning[[2]],
    #' @field alpha mixing parameter of the sparse group-penalty
    alpha = function(value) private$tuning[[3]],
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) {
      sparse <- ifelse(self$is_group_sparse, "Sparse", "")
      ridge  <- ifelse(self$is_l2_regularized, "plus L2-regularization", "")
      group  <- switch(private$type_,
             "l2" = "L1/L2 Group penalty",
             "linf" = "L1/Linf Group penalty",
             "coop" = "Cooperative Penalty"
      )
      paste(sparse, group, ridge)
    },
    #' @field group vector of integers indicating group belonging
    group = function(value) private$group_,
    #' @field type string the type of group-wise regularization applied
    type = function(value) private$type_,
    #' @field mixture_tuning mixture coefficient of the sparse group penalty
    mixture_tuning = function(value) {
      if (missing(value))
        return(private$tuning$alpha)
      else private$tuning$alpha <- value
    },
    #' @field is_group_sparse boolean indicating if sparse group or group penalty is applied
    is_group_sparse = function(value) {
      ifelse(self$mixture_tuning > 0, TRUE, FALSE)
    }
  ),
  public  = list(
    #' @description Initialize a [`SparseGroupFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param group vector of integers indicating group belonging.
    #' @param type string indicating whether the \eqn{\ell_1/\ell_2}{l1/l2} or the
    #' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted.
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, group, type, regParam) {
      super$initialize(data, intercept, regParam)
      stopifnot("The groups must be provided as a sorted vector of integers" = !is.unsorted(group))
      stopifnot("The group indices must start from 1" = min(group) == 1)
      private$group_ <- group
      private$type_  <- type
      private$optimizer <- function(dataModel, intercept, regParam, control) {
        if (data$sparse_encoding) {
          out <- switch(private$type_,
                 "linf" = group_enet_l1linf_sparse_cpp(dataModel, intercept, private$group_, regParam, control),
                 "coop" = group_enet_coop_sparse_cpp(dataModel, intercept, private$group_, regParam, control),
                 group_enet_l1l2_sparse_cpp(dataModel, intercept, private$group_, regParam, control)
          )
        } else {
          out <- switch(private$type_,
                 "linf" = group_enet_l1linf_dense_cpp(dataModel, intercept, private$group_, regParam, control),
                 "coop" = group_enet_coop_dense_cpp(dataModel, intercept, private$group_, regParam, control),
                        group_enet_l1l2_dense_cpp(dataModel, intercept, private$group_, regParam, control)
          )
        }
        out
      }
    }
  )
)

#' Class "LavaFit"
#'
#' Class of object returned by the fitting function [lava()]. Inherits fields
#' and methods of [QuadrupenFit]
#'
#' @param log_scale logical; indicates if a log-scale should be used
#' when `xvar="lambda"`. Default is `TRUE`.
#' @param standardize logical; standardize the coefficients before
#' plotting (with the norm of the predictor). Default is `TRUE`.
#' @param labels vector indicating the names associated to the plotted
#' variables. When specified, a legend is drawn in order to identify
#' each variable. Only relevant when the number of predictor is
#' small. Remind that the intercept does not count. Default is
#' \code{NULL}.
#'
#' @seealso [QuadrupenFit], [lava()]
#'
#' @export
#'
LavaFit <- R6::R6Class(
  classname = "LavaFit",
  inherit = QuadrupenFit,
  active  = list(
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) "Lava",
    #' @field lambda1 vector of tuning parameters for the l1 penalty (sparse component)
    lambda1 = function(value) private$tuning[[1]],
    #' @field lambda2 vector of tuning parameters for the l2 penalty (dense component)
    lambda2 = function(value) private$tuning[[2]],
    #' @field sparse_coef sparse part of the  decomposition of the coefficients
    sparse_coef = function(value) private$sparse_coef_,
    #' @field dense_coef dense part of the  decomposition of the coefficients
    dense_coef  = function(value) private$coef_- private$sparse_coef_,
    #' @field debias logical, should we rely on the debias coefficient of the regularizer (if available) or not
    debias = function(value) {
      if (missing(value))
        return(private$debias_)
      else {
        stopifnot(is.logical(value))
        private$debias_ <- value
        if (private$debias_) {
          private$coef_ <- private$stored_fit$coef_debias
          private$sparse_coef_ <- private$stored_fit$sparse_coef_debias
          private$intercept_ <- private$stored_fit$intercept_debias
        } else {
          private$coef_ <- private$stored_fit$coef
          private$sparse_coef_ <- private$stored_fit$sparse_coef
          private$intercept_ <- private$stored_fit$intercept
        }
      }
    }
  ),
  private = list(sparse_coef_ = NA, dense_coef_ = NA),
  public  = list(
    #' @description Initialize a [`LavaFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- lava_dense_cpp
    },
    #' @description function performing the optimization
    #' @param control list controlling the optimization process
    fit = function(control) {
      super$fit(control)
      private$sparse_coef_ <- private$stored_fit$sparse_coef
    },
    #' Plot method for lava regularization path
    #'
    #' @description Produce a plot of the solution path of a [LavaFit] object.
    #'
    #' @param xvar variable to plot on the X-axis: either `"lambda"`
    #' (\eqn{\ell_1}{l1} penalty level, or
    #' \eqn{\ell_2}{l2} for ridge  and \eqn{\ell_\infty}{l_inf}) or
    #' `"fraction"` (\eqn{\ell_1}{l1}-norm
    #' of the coefficients) or `df` for estimated degrees of freedom.
    #' Default is set to `"lambda"`.
    #' @param component a character indicating the component to plot: both
    #' (sum of sparse and dense), sparse or dense. Default to both.
    #' @param title the title. Default is set to the model name followed
    #' by what is on the Y-axis.
    #'
    #' @return a \pkg{ggplot2} object .
    plot_path = function(xvar = c("lambda", "fraction", "df"), log_scale = TRUE,
                         component = "both",
                         title = paste("Lava path:", component, "component(s)"),
                         standardize = TRUE, labels = NULL) {

      stopifnot(component %in% c("both", "sparse", "dense"))
      comp <- switch(component,
        "sparse" = self$sparse_coef, "dense" = self$dense_coef, private$coef_
        )

      xvar <- match.arg(xvar)
      stopifnot("Not available when the leading vector of tuning parameters boild down to a scalar."
                = length(self$major_tuning) > 1)

      nzeros <- which(rowSums(comp) != 0)
      stopifnot("Nothing to plot: all coefficients are zero." = length(nzeros) > 0)

      coef  <- t(as.matrix(comp[nzeros, , drop = FALSE]))
      rownames(coef) <- NULL ## avoid warning message in ggplot2
      attr(coef, "nzeros") <- nzeros

      if (standardize) coef <- scale(coef, FALSE, 1/private$data_$normx[nzeros])

      xv <- switch(xvar,
                   "fraction" = apply(abs(coef),1,sum)/max(apply(abs(coef),1,sum)),
                   "df"       = private$df_,
                   self$major_tuning
      )
      attr(xv, "type") <- xvar

      d <- .plot_regpath(xv, coef, standardize, log_scale, labels)
      d <- d + ggtitle(title) + theme_bw()
      if (is.null(labels)) d <- d + theme(legend.position = "none")
      d
    }
  )
)

#' Class "GroupLavaFit"
#'
#' Class of object returned by the fitting function [group_lava()]. Inherits fields
#' and methods of [QuadrupenFit] and [LavaFit]
#'
#' @seealso [QuadrupenFit], [group_lava()]
#'
#' @export
#'
GroupLavaFit <- R6::R6Class(
  classname = "GroupLavaFit",
  inherit = LavaFit,
  private  = list(group_ = NA, type_ = NA),
  active  = list(
    #' @field lambda1 vector of tuning parameters for the l1 group penalty (sparse component)
    lambda1 = function(value) private$tuning[[1]],
    #' @field lambda2 vector of tuning parameters for the l2 penalty (dense component)
    lambda2 = function(value) private$tuning[[2]],
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) paste0("group Lava l1/", private$type_),
    #' @field group vector of integers indicating group belonging
    group = function(value) private$group_,
    #' @field type string indicating whether the \eqn{\ell_1/\ell_2}{l1/l2} or the
    #' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted.
    type = function(value) private$type_
  ),
  public  = list(
    #' @description Initialize a [`GroupLavaFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param group vector of integers indicating group belonging.
    #' @param type string indicating whether the \eqn{\ell_1/\ell_2}{l1/l2} or the
    #' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted.
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, group, type, regParam) {
      super$initialize(data, intercept, regParam)
      stopifnot("The groups must be provided as a sorted vector of integers" = !is.unsorted(group))
      stopifnot("The group indices must start from 1" = min(group) == 1)
      private$group_ <- group
      private$type_  <- type
      private$optimizer <- function(dataModel, intercept, regParam, control) {
        out <- switch(private$type_,
                      "linf" = group_lava_l1linf_dense_cpp(dataModel, intercept, private$group_, regParam, control),
                      "coop" = group_lava_coop_dense_cpp(dataModel, intercept, private$group_, regParam, control),
                      group_lava_l1l2_dense_cpp(dataModel, intercept, private$group_, regParam, control)
        )
        out
      }
    }
  )
)
