#' Class "QuadrupenFit"
#'
#' Class of object returned by any fitting function of the
#' \pkg{quadrupen} package (\code{elastic_net} or
#' \code{bounded_reg}).
#' 
#' This class comes with the usual [predict()], [fitted()], [coef()],
#' [residuals()], [show()], [print()] and [deviance()] S3 methods.
#'
#' Specific R6 methods are available for model extraction [`QuadrupenFit$get_model()`][QuadrupenFit], 
#' cross validation [`QuadrupenFit$cross_validate()`][QuadrupenFit], stability selection 
#' [`QuadrupenFit$stability_path()`][QuadrupenFit], criteria derivation [QuadrupenFit$criteria()][QuadrupenFit] 
#' and plotting [`QuadrupenFit$plot()`][QuadrupenFit]. They come with equivalent S3 methods : [cross_validate()], 
#' [stability()] and [plot()].
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
#' @seealso See also [`InformationCriteria`], [`CrossValidation`] and
#' [`StabilityPath`]
#' 
#' @importFrom stats fitted predict residuals deviance
#'
#' @export
QuadrupenFit <- R6::R6Class(
  classname = "QuadrupenFit",
  ## ____________________________________________________
  ## 
  ## PRIVATE MEMBERS
  ## ____________________________________________________
  private = list(
    data_          = NA,
    coef_          = Matrix()  ,
    intercept_     = numeric() ,
    has_intercept_ = NA        ,
    debias_        = NA        ,
    df_            = numeric() ,
    stored_fit     = list()    ,
    tuning         = numeric() ,
    control        = list()    ,
    optimizer      = NA        ,
    infocrit       = NA        ,
    crossval       = NA        ,
    stabsel        = NA        ,
    monitoring     = list()
  ),
  ## ____________________________________________________
  ## 
  ## ACTIVE BINDINGS
  ## ____________________________________________________
  active = list(
    #' @field nvar number of coefficient (without intercept)
    nvar = function(value) {private$data_$d},
    #' @field nobs sample size
    nobs = function(value) {private$data_$n},
    #' @field dataModel an object with class [`DataModel`] storing the data
    dataModel = function(value) {private$data_},
    #' @field major_tuning vector of "leading" tuning parameters (either l1, linf or l2)
    major_tuning = function(value) {
      if (missing(value))
        return(private$tuning[[1]])
      else private$tuning[[1]] <- value
    },
    #' @field minor_tuning vector of "minor" tuning parameters (either l1 or l2)
    minor_tuning = function(value) {
      if (missing(value))
        return(private$tuning[[2]])
      else private$tuning[[2]] <- value
    },
    #' @field is_l2_regularized Boolean indicating if l2 regularization is applied
    is_l2_regularized = function(value) {
      ifelse(private$tuning$gamma > 0, TRUE, FALSE)
    },
    #' @field optim_monitoring list monitoring the optimization
    optim_monitoring = function(value) {
      if (!is.null(private$monitoring$convergence))
        private$monitoring$convergence <- 
          sapply(private$monitoring$convergence, status_to_message)
      private$monitoring
    },
    #' @field optim_config list with low level options used for optimization.
    optim_config = function(value) {private$control},
    #' @field fitted Matrix of fitted values, each column corresponding to a value of \code{lambda1}.
    fitted = function(value) {
      res <- sweep(private$data_$X %*% private$coef_ ,2L,-private$intercept_)
      res
    },
    #' @field coefficients Matrix (class `"dgCMatrix"`) of
    #' coefficients with respect to the original input. The number of
    #' rows corresponds the length of \code{lambda1}.
    coefficients         = function(value) {
      dimnames(private$coef_) <- 
        list(colnames(private$data_$X), round(c(private$tuning[[1]]),3))
      private$coef_
    },
    #' @field intercept A vector containing the successive values of the 
    #' (unpenalized) intercept.
    #' Equals to zero if \code{intercept} has been set to `FALSE`.
    intercept = function(value) {
      private$intercept_
    },
    #' @field debias logical, should we rely on the debias coefficient of the regularizer (if available) or not
    debias = function(value) {
      if (missing(value))
        return(private$debias_)
      else {
        stopifnot(is.logical(value))
        private$debias_ <- value
        if (private$debias_) {
          private$coef_ <- private$stored_fit$coef_debias
          private$intercept_ <- private$stored_fit$intercept_debias
        } else {
          private$coef_ <- private$stored_fit$coef
          private$intercept_ <- private$stored_fit$intercept
        }
      }
    },
    #' @field residuals Matrix of residuals, each column corresponding to a value of `lambda1`.
    residuals       = function(value) {apply(self$fitted, 2, function(y_hat) private$data_$y - y_hat)},
    #' @field deviance the model deviance
    deviance        = function(value) {colSums(self$residuals^2)},
    #' @field degrees_freedom Estimated degree of freedoms for the successive `lambda1`.
    degrees_freedom = function(value) {private$df_},
    #' @field r_squared vector giving the coefficient of determination as a function of lambda1.
    r_squared       = function(value) {
      rss <- sum((private$data_$y - ifelse(private$has_intercept_, mean(private$data_$y), 0))^2)
      1 - colSums(self$residuals^2) / rss
    },
    #' @field information_criteria object with class [`InformationCriteria`] storing various information criteria 
    #' (AIC, BIC, GCV, etc) for the current fit.
    information_criteria = function(value) {private$infocrit},
    #' @field cross_validation object with class [`CrossValidation`] storing output of CV job. 
    #' Only available once method cross_validate has been called.
    cross_validation = function(value) {private$crossval},
    #' @field stability_path object with class [`StabilityPath`] storing output of stability selection. 
    #' Only available once method $stability has been called.
    stability_path   = function(value) {private$stabsel}
  ),
  
  ## ____________________________________________________
  ## 
  ## PUBLIC MEMBERS
  ## ____________________________________________________
  public  = list(
    
    #' @description Initialize a [QuadrupenFit] model
    #' @param data a [DataModel] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize = function(data, intercept, regParam) {

      stopifnot("The data object must be an instance of DataModel" = inherits(data, "DataModel"))
      stopifnot("regParam must be a list" = is.list(regParam))
      stopifnot("All regularization parameters must be positive." = all(unlist(regParam) >= 0))
      stopifnot("The first entry of regParam must be sorted in decreasing order." =
                  !is.unsorted(rev(regParam[[1]])))
      stopifnot("The second entry of regParam must be a scalar (cannot be a vector)." = 
                  (length(regParam[[2]]) == 1 & inherits(regParam[[2]], "numeric")))
      stopifnot("minratio must be non negative." = regParam$min_ratio > 0)
      stopifnot("nlambda1 must be non negative." = regParam$n_lambda1 > 0)
      private$data_          <- data
      private$has_intercept_ <- intercept
      private$tuning         <- regParam
      private$debias_        <- FALSE
    },
    #' @description User friendly print method
    show = function() {
      cat("Linear regression with", self$penalty, "penalizer.\n")
      cat("- number of coefficients:", self$nvar,"+ intercept\n")
      cat("- ", names(private$tuning)[[1]], " regularization: ",
          length(self$major_tuning), " points from ",
          format(max(self$major_tuning), digits = 3)," to ",
          format(min(self$major_tuning), digits = 3),"\n", 
          "- ", names(private$tuning)[[2]], " regularization: ",
          self$minor_tuning, "\n", sep = ""
        )
      invisible(self)
    },
    #' @description User friendly print method
    print = function() { self$show() },
    #' @description function performing the optimization
    #' @param control list controlling the optimization process
    fit = function(control) {
      ## ======================================================
      ## C++ CALL OPTIMIZER
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      private$stored_fit <- private$optimizer(private$data_, private$has_intercept_, private$tuning, control)
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA)
      ## END OF CALL
      ## ======================================================
      private$tuning[[1]] <- private$stored_fit$tuning_param[[1]]
      private$intercept_  <- drop(private$stored_fit$intercept)
      private$coef_       <- private$stored_fit$coef
      private$df_         <- drop(private$stored_fit$df)
      private$monitoring  <- private$stored_fit$monitoring
      private$monitoring$timer <- timer
      private$control     <- control
    },
    #' @description Model extraction
    #' @param selection either a character (model selection criteria) of a scalar (lambda value)
    #' @param type character for the desired output
    #' @return either a vector of coefficients, a scalar or the model index
    get_model = function(
          selection,
          type = c("coefficients", "penalty", "index")) {
      lambda <- private$tuning[[1]]
      if (is.character(selection)) {
        stopifnot("must be a character in" = selection %in% c("AIC", "BIC", "mBIC", "eBIC", "GCV", "CV_min", "CV_1se"))
        if (selection %in% c("CV_min", "CV_1se") & !inherits(private$crossval, "CrossValidation")) 
          stop("Cross-validation has not yet been performed")
        
        index <- 
          switch(selection,
                 "AIC"  = which.min(private$infocrit$data$AIC),
                 "BIC"  = which.min(private$infocrit$data$BIC),
                 "mBIC" = which.min(private$infocrit$data$mBIC),
                 "eBIC" = which.min(private$infocrit$data$eBIC),
                 "GCV"  = which.min(private$infocrit$data$GCV),
                 "CV_min" = min(match(private$crossval$lambda1_min, lambda), length(lambda)),
                 "CV_1se" = min(match(private$crossval$lambda1_1se, lambda), length(lambda)),
          )
      } else {
        index <- match(selection, lambda)
        if (is.na(index)) { ## No exact match
          index <- which.min(abs(selection - lambda)) ## closest model (in terms of parameter value)
          warning(paste("No such a model in the collection. Acceptable parameter values can be found via $major_tuning",
                        paste0("  Returning model with closest value. Requested: ", selection, ", returned: ", lambda[index]),
                        sep = "\n"))
        }
      }
      type <- match.arg(type)
      if (type == "index")
        return(index)
      else if (type == "penalty")
        return(lambda[index])
      else {
        if (private$has_intercept_) {
          res <- setNames(c(private$intercept_[index], private$coef_[,index]), c("Intercept", private$data_$varnames))
        } else {
          res <- setNames(private$coef_[,index], private$data_$varnames)
        }
        return(res)
      }
    },
    #' @description Predict response for new sample based on the current model
    #' @param newx matrix of new values for the regressor with which to predict. If omitted, the fitted values are used.
    #' @param selection either a character (model selection criteria) of a scalar (lambda value)
    #' @return a vector of predicted value
    predict = function(newx = NULL, selection = NULL) {
      
      if (is.null(selection)) {
        index <- 1:length(private$tuning[[1]])
      } else {
        index <- self$get_model(selection, type = "index")
      }

      if (is.null(newx)) {
        res <- self$fitted[ , index, drop = FALSE]
      } else {
        res <- sweep(newx %*% private$coef_[, index , drop = FALSE], 2L, -private$intercept_[index])
      }
      res
    },
    #' Cross-validation for Quadrupen object
    #' 
    #' @description Function that computes K-fold cross-validated error of a
    #' \code{quadrupen} fit, possibly on a grid of `lambda1`, `lambda2`.
    #'
    #' @param K integer indicating the number of folds. Default is 10.
    #'
    #' @param folds list of `K` vectors that describes the folds to
    #' use for the cross-validation. By default, the folds are randomly
    #' sampled with the specified K. The same folds are used for each
    #' values of `lambda2`.
    #'
    #' @param lambda2 tunes the \eqn{\ell_2}{l2}-penalty (ridge-like) of
    #' the fit. If none is provided, a vector of values is generated and
    #' a CV is performed on a grid of `lambda2` and `lambda1`,
    #' using the same folds for each `lambda2`.
    #'
    #' @param verbose logical; indicates if the progression (the current
    #' `lambda2` should be displayed. Default is `TRUE.`
    #'
    #' @param cores the number of cores to use. The default uses 1 core 
    #' (safer in case your BLAS/LAPACK libraries are multithreaded)
    #'
    #' @return an object with class [CrossValidation] is sent back and stored as a 
    #' field of the original [QuadrupenFit] object.
    #' 
    #' @seealso [cross_validate()]
    cross_validate = 
      function(
          K       = 10,
          folds   = split(sample(1:self$nobs), rep(1:K, length = self$nobs)),
          lambda2 = self$minor_tuning, verbose = TRUE, cores = 1) {

        ## Some variables and copies useful for CV work
        K <- length(folds)
        control  <- private$control; control$verbose <- 0 ; control$adjust <- FALSE
        regParam <- private$tuning # copy of the currently used tuning parameters
        nlambda1 <- length(regParam[[1]])
        nlambda2 <- length(lambda2)
        lambda2_vec <- rep(lambda2, each = K)
        fold_id <- rep(1:K, nlambda2)
        
        if (verbose) {
          cat("\nCROSS-VALIDATION FOR ", self$penalty," REGULARIZER \n\n")
          cat(K, "-fold CV on a grid of (", 
              nlambda1, ",", nlambda2, ") tuning parameters\n", sep = "")
        }

        ## Same data splitting is kept for varying lambda2 values
        CVData <- self$dataModel$splitTrainTest(K, folds)

        ## CV err for a fixed couple fold/lambda2
        one_fold <- function(fold, lambda2) {
          if (verbose & (fold == 1)) cat(round(lambda2, 3),"\t")
          regParam[[2]] <- lambda2
          out <- private$optimizer(CVData[[fold]]$trainData, private$has_intercept_, regParam, control)
          if (private$debias_) {
            intercept   <- out$intercept_debiased
            coef <- out$coef_debiased
          } else {
            intercept   <- out$intercept
            coef <- out$coef
          }
          y_hat <- scale(CVData[[fold]]$testData$X %*% coef, -intercept, FALSE)
          err <- sweep(y_hat, 1L, CVData[[fold]]$testData$y)^2
          if (ncol(err) < length(regParam[[1]])) {
            NAs <- length(regParam[[1]]) - ncol(err)
            err <- cbind(err, matrix(NA, nrow(err),NAs))
          }
          err
        }

        err <- do.call(rbind, 
          parallel::mcmapply(FUN = one_fold, fold = fold_id, lambda2 = lambda2_vec, 
                   mc.cores = cores,
                   mc.preschedule = ifelse(K > 10,TRUE,FALSE),
                   SIMPLIFY = FALSE
          )) |> as.matrix() |> as.data.frame()
        if (verbose) cat("\n")

        res <- do.call(rbind, tapply(err, rep(1:nlambda2, each = self$nobs), function(err_) {
          mean <- colMeans(err_, na.rm = TRUE)
          if (any(is.nan(mean))) {
            warning("\nThere have been a lot of early stops along the path: 
                  I keep on running, but you should reconsider the value of 
                  the minimal penalty along the path regarding the n<<p setting.")
          }
          mean[is.nan(mean)] <- NA
          serr <- colSums(sweep(err_, 2L, mean, check.margin = FALSE)^2, na.rm = TRUE)
          serr <- sqrt((serr/(self$dataModel$n - 1))/K)
          data.frame(mean = mean, serr = serr, lambda1 = private$tuning[[1]])
        }))
        res$lambda2 <- rep(lambda2, each = length(private$tuning[[1]]))

        private$crossval <- CrossValidation$new(cv_error = res, folds = folds)
        invisible(private$crossval)
    },
    #' Stability selection for Quadrupen object
    #' 
    #' @description Compute the stability path of a (possibly randomized) fitting
    #' procedure as introduced by Meinshausen and Buhlmann (2010).
    #' 
    #' @param n_subsamples integer indicating the number of subsamplings
    #' used to estimate the selection probabilities. Default is 100.
    #'
    #' @param subsample_size integer indicating the size of each subsamples.
    #' Default is `floor(n/2)`.
    #'
    #' @param subsamples list with `subsamples` entries with vectors
    #' describing the folds to use for the stability procedure. By
    #' default, the folds are randomly sampled with the specified
    #' \code{n_subsamples} and `subsample_size` argument.
    #'
    #' @param weakness Coefficient used for randomizing the weights of each features.
    #' Default is 1` for no randomization. See details below.
    #' 
    #' @param verbose logical; indicates if the progression should be
    #' displayed. Default is `TRUE`.
    #'
    #' @param cores the number of cores to use. The default uses 1 core 
    #' (safer in case your BLAS/LAPACK libraries are multithreaded)
    #'
    #' @return an object with class [StabilityPath] is sent back and stored as a 
    #' field of the original [QuadrupenFit] object.
    #' 
    #' @seealso [stability()]
    stability = function(
          n_subsamples   = 50,
          subsample_size = floor(self$nobs/2),
          subsamples     = replicate(n_subsamples, sample(1:self$nobs, subsample_size), simplify = FALSE),
          weakness       = 1,
          verbose        = TRUE,
          cores          = 1) {
      
      ## =============================================================
      ## INITIALIZATION & PARAMETERS RECOVERY
      if (Sys.info()[['sysname']] == "Windows") {
        warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
        cores <- 1
      }
      
      control  <- private$control; 
      control$verbose <- 0 ; control$adjust <- FALSE
      control$maxfeat  <- self$nvar
      nlambda1 <- length(private$tuning[[1]])

      ## Prepare blocs of sub samples to run jobs parallely
      blocs <- suppressWarnings(split(1:n_subsamples, 1:cores))
      
      if (verbose) {
        cat(paste("\n\nSTABILITY SELECTION ",ifelse(weakness < 1,"with","without")," randomization (weakness = ",weakness,")",sep = ""))
        cat(paste("\nFitting procedure: ", self$penalty," with ", nlambda1,"-dimensional grid of lambda1.", sep = ""))
        cat("\nRunning",length(blocs),"jobs parallely (1 per core)")
        cat("\nApprox.", length(blocs[[1]]),"subsamplings for each job for a total of", n_subsamples)
      }
      ## get data samples
      SubsampledData <- self$dataModel$splitSubSamples(n_subsamples, subsample_size, subsamples, weakness)

      ## function to run on each core
      bloc_stability <- function(subsets) {
        select <- Matrix(0, self$nvar, nlambda1)
        subsamples_ok <- 0
        for (s in 1:length(subsets)) {
          active <- private$optimizer(SubsampledData[[subsets[s]]], private$has_intercept_, private$tuning, control)$active
          if (ncol(active) == nlambda1) {
            subsamples_ok <- subsamples_ok + 1
            select <- select + active
          }
        }
        if (subsamples_ok < 0.5*length(subsets)) {
          cat("\nWarning: more than 50% of the run were discarded in that core due to early stops of the fitting procedure. You should consider largest 'minratio' or strongest 'lambda2'.")
        }
        return(select/(subsamples_ok*length(blocs)))
      }
      
      ## Now launch the B jobs...
      prob_bloc <- mclapply(blocs, bloc_stability, mc.cores = cores)
      
      ## Construct the probability path
      path <- Matrix(0, self$nvar, nlambda1)
      for (b in 1:length(prob_bloc)) {
        path <- path + prob_bloc[[b]]
      }

      private$stabsel <- StabilityPath$new(
        probabilities = path          ,
        regParam      = private$tuning,
        subsamples    = subsamples
      )
      invisible(private$stabsel)
    },
    #' Penalized criteria based on estimation of degrees of freedom
    #'
    #' @description Produce a plot or send back the values of some penalized criteria
    #' accompanied with the vector(s) of parameters selected
    #' accordingly. The default behavior plots the BIC and the AIC (with
    #' respective factor \eqn{\log(n)}{log(n)} and \eqn{2}{2}) yet the user can specify any
    #' penalty.
    #'
    #' @param penalty a vector with as many penalties a desired. The
    #' default contains the penalty corresponding to the AIC and the BIC
    #' (\eqn{2}{2} and \eqn{\log(n)}{log(n)}). Setting the "names"
    #' attribute, as done in the default definition, leads to outputs
    #' which are easier to read.
    #' @param sigma scalar: an estimate of the residual variance. When
    #' available, it is plugged-in the criteria, which may be more
    #' relevant. If `NULL` (the default), it is estimated as usual
    #' (see details).
    #' 
    #' @return an object with class [InformationCriteria] is sent back and stored as a 
    #' field of the original [QuadrupenFit] object.
    #' 
    #' @seealso [criteria()]
    criteria = function(penalty=
                          setNames(c(2, log(self$nobs), log(self$nvar), log(self$nobs) + 2*log(self$nvar)),
                                   c("AIC","BIC", "mBIC", "eBIC")), sigma=NULL) {
      
      if (is.null(sigma)) {
        crit <- sapply(penalty, function(pen) self$nobs*log(self$deviance/self$nobs) + pen * self$degrees_freedom)
      } else {
        crit <- sapply(penalty, function(pen) self$deviance/sigma^2 + pen * self$degrees_freedom)
      }
      crit <- as.data.frame(crit)
      ## Compute generalized cross-validation
      crit$GCV <- self$deviance/(self$nobs*(1 - self$degrees_freedom/self$nobs)^2)
      
      ## Put together all relevant information about those criteria
      private$infocrit <- 
        InformationCriteria$new(
          value = data.frame(
            crit, 
            df        = self$degrees_freedom, 
            lambda    = self$major_tuning, 
            fraction  = colSums(abs(private$coef_))/max(colSums(abs(private$coef_))), 
            row.names = 1:nrow(crit)
          )
        )
      invisible(private$infocrit)
    },
    #' @description Plot method for QuadrupenFit
    #' 
    #' @param type the type of plot, either `"path"` for regularization path; 
    #' `"criteria"` for BIC-like  information criteria ; `"crossval"` for 
    #' cross-validation plot ; and `"stability"` for stability path.
    #' 
    #' @details  The `"path"`plot is available as soon as a fit has been performed.
    #' For the others, the appropriate post-treatments must have been made via the
    #' methods [`QuadrupenFit$criteria()`][QuadrupenFit], [`QuadrupenFit$cross_validate()`][QuadrupenFit] or
    #' [`QuadrupenFit$stability()`][QuadrupenFit]
    #' 
    #' All plots functions are given with the default arguments, except for `labels` and `log_scale`.
    #' If you need more control, please use the dedicated methods: [`QuadrupenFit$plot_path()`][QuadrupenFit], 
    #' [`InformationCriteria$plot()`][InformationCriteria], [`CrossValidation$plot()`][CrossValidation],
    #'  [`StabilityPath$plot()`][StabilityPath] or the corresponding S3 methods.
    #' 
    plot = function(type = c("path", "criteria", "crossval", "stability"), 
                    log_scale = TRUE, labels = NULL) {
      
      type <- match.arg(type)

      stopifnot("Not available: the method $criteria() has not been called yet." = 
        (inherits(self$information_criteria, "InformationCriteria") | type != "criteria"))
      stopifnot("Not available: the method $cross_validate() has not been called yet." = 
        (inherits(self$cross_validation, "CrossValidation") | type != "crossval"))
      stopifnot("Not available: the method $stability() has not been called yet." = 
        (inherits(self$stability_path, "StabilityPath") | type != "stability"))
      
      d <- switch(type, 
              "path"      = self$plot_path(log_scale = log_scale, labels = labels),
              "criteria"  = self$information_criteria$plot(log_scale = log_scale),
              "crossval"  = self$cross_validation$plot(log_scale = log_scale),
              "stability" = self$stability_path$plot(labels = labels)
      )
      d
      
    },
    #' Plot method for regularization path
    #' 
    #' @description Produce a plot of the solution path of a [QuadrupenFit] object.
    #'
    #' @param xvar variable to plot on the X-axis: either `"lambda"`
    #' (\eqn{\ell_1}{l1} penalty level, or
    #' \eqn{\ell_2}{l2} for ridge  and \eqn{\ell_\infty}{l_inf}) or
    #' `"fraction"` (\eqn{\ell_1}{l1}-norm
    #' of the coefficients) or `df` for estimated degrees of freedom. 
    #' Default is set to `"lambda"`.
    #' @param title the title. Default is set to the model name followed
    #' by what is on the Y-axis.
    #'
    #' @return a \pkg{ggplot2} object .
    #'
    #' @examples \dontrun{
    #' ## Simulating multivariate Gaussian with blockwise correlation
    #' ## and piecewise constant vector of parameters
    #' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
    #' cor <- 0.75
    #' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
    #' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
    #' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
    #' diag(Sigma) <- 1
    #' n <- 50
    #' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
    #' y <- 10 + x %*% beta + rnorm(n,0,10)
    #'
    #' ## Plot the Lasso path
    #' plot(lasso(x,y), title="Lasso solution path")
    #' ## Plot the Elastic-net path
    #' plot(elastic_net(x,y), title = "Elastic-net solution path")
    #' ## Plot the Elastic-net path (fraction on X-axis, unstandardized coefficient)
    #' plot(elastic_net(x,y, lambda2=10), standardize=FALSE, xvar="fraction")
    #' ## Plot the Bounded regression path (fraction on X-axis)
    #' plot(bounded_reg(x,y, lambda2=10), xvar="fraction")
    #' }
    #'
    plot_path = function(xvar = c("lambda", "fraction", "df"), log_scale = TRUE,
                    title = paste("Path for", self$penalty),
                    standardize = TRUE, labels = NULL) {

      xvar <- match.arg(xvar)      
      stopifnot("Not available when the leading vector of tuning parameters boild down to a scalar." 
                = length(self$major_tuning) > 1)

      nzeros <- which(rowSums(private$coef_) != 0)
      stopifnot("Nothing to plot: all coefficients are zero." = length(nzeros) > 0)
      
      coef  <- t(as.matrix(private$coef_[nzeros, , drop = FALSE]))
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
