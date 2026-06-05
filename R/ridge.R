#' Fit a linear model with a structured ridge regularization
#'
#' Adjust a linear model with ridge regularization (possibly
#' structured \eqn{\ell_2}{l2}-norm). The solution path is computed
#' at a grid of values for the \eqn{\ell_2}{l2}-penalty. See details
#' for the criterion optimized.
#'
#' @inheritParams elastic_net
#' 
#' @param lambda sequence of decreasing \eqn{\ell_2}{l2}-penalty
#' levels. If `NULL` (the default), a vector is generated with
#' `nlambda` entries, starting from a guessed level
#' `lambda_max` where only the intercept is included, then
#' shrunken to `minratio*lambda_max`.
#' 
#' @param nlambda integer that indicates the number of values to put
#' in the `lambda` vector.  Ignored if `lambda` is provided.
#'
#' @param lambda_max the largest value of `lambda` considered
#' 
#' @return an object with class [RidgeRegressionFit], inheriting from [QuadrupenFit].
#'
#' @details The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda_2} = \arg \min_{\beta} \frac{1}{2} (y - X
#' \beta)^T (y - X \beta) + \frac{\lambda_2}{2} \beta^T S \beta, }}
#' \if{html}{\out{  &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>2</sub></sub> = argmin<sub>&beta;</sub> 1/2
#' RSS(&beta;) + &lambda;/2 <sub>2</sub> &beta;<sup>T</sup> S
#' &beta;,  }} \if{text}{\deqn{beta.hat(lambda2) =
#' argmin_beta 1/2 RSS(beta) + lambda2 beta' S beta,}} where the
#' \eqn{\ell_2}{l2} structuring positive semidefinite matrix
#' \eqn{S}{S} is provided via the \code{struct} argument (possibly of
#' class \code{Matrix}).
#'
#' @seealso See also [QuadrupenFit]
#' 
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' cor <- 0.75
#' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
#' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
#' Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo)
#' diag(Sigma) <- 1
#' n <- 50
#' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
#' y <- 10 + x %*% beta + rnorm(n,0,10)
#'
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- "relevant"
#' plot(ridge(x,y) , label=labels) ## a mess
#' plot(ridge(x,y, struct=solve(Sigma)), label=labels) ## even better
#'
#' @export
ridge <- function(x,
                  y,
                  lambda     = NULL,
                  weights    = rep(1, nrow(x)),
                  struct     = Matrix::Diagonal(ncol(x),1),
                  penscale   = rep(1,ncol(x)),
                  intercept  = TRUE,
                  normalize  = TRUE,
                  nlambda    = 100 ,
                  minratio   = 1e-5,
                  lambda_max = 100,
                  control    = list()) {

  ## ============================================
  ## RECOVER LOW LEVEL OPTIONS
  ## 
  ctrl <- list(verbose = 1, # default control options
               timer   =  FALSE)
  ctrl[names(control)] <- control # overwritten by user specifications
  ctrl$normalize <- normalize
  
  ## ============================================
  ## INSTANTIATE THE DATA MODEL
  ## 
  myData <- DataModel$new(
    covariates  = as.matrix(x),
    outcome     = y,
    cov_struct  = struct,
    obs_weights = weights
  )
  myData$CholStruct()

  ## ============================================
  ## INSTANTIATE THE PENALTY MODEL
  ## 
  myModel <- RidgeRegressionFit$new(
    data      = myData,
    intercept = intercept,
    regParam  = list(lambda = lambda, 
                     gamma  = 0,
                     lambda_factor = penscale, 
                     min_ratio = minratio,
                     n_lambda = nlambda)
  )

  ## ============================================
  ## FIT THE MODEL WITH ACTIVE SET ALGORITHM
  myModel$fit(ctrl)

  ## ============================================
  ## POSTREATMENT + SEND BACK THE RESULTING MODEL
  ##
  myModel$criteria()
  myModel
}
