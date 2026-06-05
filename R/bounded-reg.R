#' Fit a linear model with infinity-norm plus ridge-like regularization
#'
#' Adjust a linear model penalized by a (possibly
#' weighted) \eqn{\ell_\infty}{l-infinity}-norm (bounding the
#' magnitude of the parameters) and a (possibly structured)
#' \eqn{\ell_2}{l2}-norm (ridge-like). The solution path is computed
#' at a grid of values for the infinity-penalty, fixing the amount of
#' \eqn{\ell_2}{l2} regularization. See details for the criterion
#' optimized.
#'
#' @inheritParams sparse_lm
#' 
#' @return an object with class [QuadrupenFit].
#'
#' @details The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda_1,\lambda_1} = \arg \min_{\beta} \frac{1}{2}
#' (y - X \beta)^T (y - X \beta) + \lambda_1 \|D \beta \|_{\infty} +
#' \frac{\lambda_2}{2} \beta^T S \beta, }} \if{html}{\out{ 
#' &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
#' argmin<sub>&beta;</sub> 1/2 RSS(&beta;) + &lambda;<sub>1</sub>
#' &#124; D &beta; &#124;<sub>&infin;</sub> + &lambda;/2 <sub>2</sub>
#' &beta;<sup>T</sup> S &beta;,  }}
#' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
#' RSS(beta) + lambda1 max|D beta| + lambda2 beta' S beta,}} where
#' \eqn{D}{D} is a diagonal matrix, whose diagonal terms are provided
#' as a vector by the \code{penscale} argument. The \eqn{\ell_2}{l2}
#' structuring matrix \eqn{S}{S} is provided via the \code{struct}
#' argument, a positive semidefinite matrix (possibly of class
#' \code{Matrix}).
#'
#' Note that the quadratic algorithm for the bounded regression may
#' become unstable along the path because of singularity of the
#' underlying problem, e.g. when there are too much correlation or
#' when the size of the problem is close to or smaller than the
#' sample size. In such cases, it might be a good idea to switch to
#' the proximal solver, slower yet more robust. This is the strategy
#' automatically adopted in code, that will send a warning in verbose mode
#' while switching the method to \code{'fista'} and keep on
#' optimizing on the remainder of the path.
#'
#' Singularity of the system can also be avoided with a larger
#' \eqn{\ell_2}{l2}-regularization, via \code{lambda2}, or a
#' "not-too-small" \eqn{\ell_\infty}{l-infinity} regularization, via
#' a larger \code{'minratio'} argument.
#'
#' @return an object with class [BoundedRegressionFit], inheriting from [QuadrupenFit].
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
#' ## Infinity norm without/with an additional l2 regularization term
#' ## and with structuring prior
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- "relevant"
#' plot(bounded_reg(x,y,lambda2=0) , label=labels) ## a mess
#' plot(bounded_reg(x,y,lambda2=10), label=labels) ## good guys are at the boundaries
#' plot(bounded_reg(x,y,lambda2=10,struct=solve(Sigma)), label=labels) ## even better
#'
#' @export
bounded_reg <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.01,
                        penscale  = rep(1,ncol(x)),
                        struct    = Matrix::Diagonal(ncol(x), 1),
                        intercept = TRUE,
                        normalize = TRUE,
                        nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                        minratio  = ifelse(nrow(x) <= ncol(x), 1e-2, 1e-4),
                        maxfeat   = ifelse(lambda2 < 1e-2, min(nrow(x),ncol(x)), min(4*nrow(x),ncol(x))),
                        control   = list()) {
  
  ## ============================================
  ## RECOVER LOW LEVEL OPTIONS
  ## 
  ctrl <- optim_breg_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  # if (!is.null(control$method)) if (control$method != "quadra") ctrl$threshold <- 1e-2
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$method <- switch(ctrl$method, quadra = "QUADRA", pgd = "PGD", fista = "FISTA", 0)
  ctrl$normalize <- normalize
  
  ## ============================================
  ## INSTANTIATE THE DATA MODEL
  ## 
  myData <- DataModel$new(
    covariates  = x,
    outcome     = y,
    cov_struct  = struct
  )
  
  ## ============================================
  ## INSTANTIATE THE PENALTY MODEL
  myModel <- BoundedRegressionFit$new(
    data      = myData,
    intercept = intercept,
    regParam  = list(lambda = lambda1,
                     gamma  = lambda2, 
                     lambda_factor = penscale,
                     min_ratio = minratio, n_lambda = nlambda1)
  )
  
  ## ============================================
  ## FIT THE MODEL WITH ACTIVE SET ALGORITHM
  ##
  if (ctrl$verbose) cat("\nModel fitting and optimization")
  myModel$fit(ctrl)

  ## ============================================
  ## POSTREATMENT + SEND BACK THE RESULTING MODEL
  ##
  if (ctrl$verbose) cat("\nPost-treatment")
  myModel$criteria()
  myModel
}

#' @rdname bounded_reg
#' @importFrom lifecycle badge deprecate_warn
#' @export
bounded.reg <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.01,
                        penscale  = rep(1,ncol(x)),
                        struct    = Matrix::Diagonal(ncol(x), 1),
                        intercept = TRUE,
                        normalize = TRUE,
                        nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                        minratio  = ifelse(nrow(x) <= ncol(x), 1e-2, 1e-4),
                        maxfeat   = ifelse(lambda2 < 1e-2, min(nrow(x),ncol(x)), min(4*nrow(x),ncol(x))),
                        control   = list()) {

  lifecycle::deprecate_warn("1.1.0", "bounded.reg()", "bounded_reg()")
  
  out <- bounded_reg(x,
                     y,
                     lambda1   = lambda1,
                     lambda2   = lambda2,
                     penscale  = penscale,
                     struct    = struct,
                     intercept = intercept,
                     normalize = normalize,
                     nlambda1  = nlambda1,
                     minratio  = minratio,
                     maxfeat   = maxfeat,
                     control   = control)
  out
}

                        
