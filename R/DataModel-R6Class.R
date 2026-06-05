#' Data Class 
#' 
#' @description Class for storing data and various fixed quantity
#' 
#' @export
DataModel <- R6::R6Class(
  classname = "DataModel",
  private = list(
    names     = NA
  ),
  public = list(
    ## model-related fields
    #' @field X matrix of regressor 
    X = Matrix(),
    #' @field y vector of response
    y = numeric(),
    #' @field C_inv Inverse of the Cholesky decomposition of S
    C_inv = matrix(),
    #' @field S SDP structuring matrix 
    S = Matrix(),
    #' @field wy vector of observation weights
    wy = numeric(),
    initialize = 
      #' @description constructor for DataModel
      #' @param covariates matrix of covariates/regressors
      #' @param outcome vector of outcome/response
      #' @param cov_struct sdp matrix structuring the covariates/regressors
      #' @param cov_weights vector of covariates/regressors weights
      #' @param obs_weights vector of observations weights
      #' @param check_args logical, should args be check at initialization?
      function(covariates, outcome, cov_struct,
               obs_weights = rep(1,length(outcome)), check_args = TRUE) {
        
        ## ===================================================
        ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
        ##
        if (check_args) {
          stopifnot("x has to be of class 'matrix' or 'dgCMatrix'." = 
                      inherits(covariates, c("matrix", "dgCMatrix")))
          stopifnot("NA value in x not allowed." = !any(is.na(covariates)))
          stopifnot("y has to be of type 'numeric'" = is.numeric(outcome))
          stopifnot("x and y have not correct dimensions" = 
                      (nrow(covariates) == length(outcome)))
          stopifnot("cov_struct does not have the appropriate type" = 
            any(c(inherits(cov_struct, "matrix"), inherits(cov_struct, "sparseMatrix"), 
                  inherits(cov_struct, "FusionGraph"))))
          if (is.matrix(cov_struct) | inherits(cov_struct, "sparseMatrix")) {
            cov_struct <- as(cov_struct, "CsparseMatrix")
            stopifnot("struct must be a (square) positive semidefinite matrix." = 
                        all(dim(cov_struct) == ncol(covariates)))
          }
          stopifnot("observations weights must be positive" = all(obs_weights > 0))
        } 
        ## ===================================================
        
        if (is.null(colnames(covariates))) colnames(covariates) <- 1:ncol(covariates)
        self$X  <- covariates
        self$y  <- outcome
        self$S  <- cov_struct
        self$wy <- obs_weights
      },
    #' @description Compute Cholesky factorization of the Structuring matrix 
    CholStruct = function() {
      stopifnot(inherits(self$S, "sparseMatrix"))
      self$C_inv <- as.matrix(solve(Matrix::chol(self$S)))
    },
    #' @description a function splitting the data into train and test folds
    #' @param nfolds the number of folds
    #' @param folds a list of vectors describing the folds (optional)
    #' @return a list with train and test data and id.
    splitTrainTest = function(
    nfolds = 10, folds  = split(sample(1:self$n), rep(1:nfolds, length = self$n))
    ) {
      ## create the list of split each compose with couple of Train/Test
      lapply(folds, function(omit) {
        trainData <- DataModel$new(self$X[-omit, , drop = FALSE], self$y[-omit], self$S, self$wy[-omit])
        testData  <- DataModel$new(self$X[ omit, , drop = FALSE], self$y[omit], self$S, self$wy[omit])
        trainData$C_inv <- self$C_inv # Cholesky factorization remain the same
        testData$C_inv  <- self$C_inv # 
        list(trainData = trainData, testData = testData,
             trainID = setdiff(1:self$n, omit), testID = omit)
      })
    },
    #' @description a function splitting data into subsamples
    #' @param n_subsamples the number of subsamples
    #' @param subsample_size the subsample size
    #' @param subsamples list with vector of subsamples (optional)
    #' @param weakness coefficient for randomly weighting the regressor, default to 1
    #' @return  a list of DataModel, resampling of the original
    splitSubSamples = function(
    n_subsamples = 50,
    subsample_size = floor(self$n/2),
    subsamples = replicate(n_subsamples, sample(1:self$n, subsample_size), simplify = FALSE),
    weakness = 1
    ) {
      ## create the list of split each compose with couple of Train/Test
      lapply(subsamples, function(keep) {
        Xs <- Matrix::colScale(self$X[keep, ], runif(self$d, weakness, 1))
        if (!inherits(self$X, "sparseMatrix")) Xs <- as.matrix(Xs)
        DataModel$new(Xs, self$y[keep], self$S, self$wy[keep])
      })
    }
  ), 
  active = list(
    #' @field d number of regressor
    d = function() ncol(self$X),
    #' @field n sample size
    n = function() nrow(self$X),
    #' @field sparse_encoding logical indicating if the matrix of regressor is sparsely encoded
    sparse_encoding = function() {inherits(self$X, "sparseMatrix")},
    #' @field varnames character, the names of the covariates/regressors
    varnames = function() {colnames(self$X)},
    #' @field normx norm of each column of X
    normx = function() {sqrt(drop(colSums(self$X^2)) - self$n * colMeans(self$X)^2)}
  )
)
