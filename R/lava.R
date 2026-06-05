#' Fit a linear model with lava regularization
#'
#' Adjust a lava regularized linear model, that is a lava transformation 
#' of the data followed by a 
#' (possibly weighted) \eqn{\ell_1}{l1}-norm. The solution path is
#' computed at a grid of values for the \eqn{\ell_1}{l1}-penalty. See
#' details for the criterion optimized.
#'
#' @inheritParams elastic_net
#'
#' @return an object with class [QuadrupenFit].
#'
#' @details The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\theta}_{\lambda_1,\lambda_2} = \arg \min_{\theta = \beta + \delta} \frac{1}{2 n} (y - X
#' (\beta + \delta))^T (y - X (\beta + \delta)) + \lambda_1 \| \beta \|_{1} + \frac{\lambda_2}{2} \delta^T S \delta, }}
#' \if{html}{\out{  &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub></sub> =
#' argmin<sub>&theta; = &beta;+&delta;</sub> 1/2n RSS(&beta; + &delta;) + &lambda;<sub>1</sub>
#' &#124; &beta; &#124;<sub>1</sub> + &lambda;/2 <sub>2</sub> &delta;<sup>T</sup> S
#' &delta;,  }}
#' \if{text}{\deqn{theta.hat(lambda1,lambda2) = argmin_{theta = beta + delta} 1/2n
#' RSS(beta + delta) + lambda1 |beta|1 + + lambda2 delta' S delta}}.
#'
#' @return an object with class [LavaFit], inheriting from [QuadrupenFit].
#' 
#' @references Chernozhukov, Victor, Christian Hansen, and Yuan Liao. "A lava attack on the recovery of sums of 
#' dense and sparse signals." The Annals of Statistics (2017): 39-76. <doi:10.1214/16-AOS1434>
#' 
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta  <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' delta <- runif(sum(c(25,10,25,10,25)),-.1,.1)
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
#' ## The solution path of the LAVA
#' out <- lava(x,y)
#' out$plot_path(component = "sparse", labels=labels)
#' out$plot_path(component = "dense", labels=labels)
#'
#' @export
lava <- function(x,
                 y,
                 lambda1   = NULL,
                 lambda2   = 1,
                 weights   = rep(1, nrow(x)),
                 penscale  = rep(1,ncol(x)),
                 struct    = Matrix::Diagonal(ncol(x), 1),
                 intercept = TRUE,
                 normalize = TRUE,
                 refit     = FALSE,
                 nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                 minratio  = 1e-2,
                 maxfeat   = min(nrow(x),ncol(x)),
                 beta0     = numeric(ncol(x)),
                 control   = list()) {
  
  ## ============================================
  ## RECOVER LOW LEVEL CONFIGURATION
  ##
  ctrl <- optim_enet_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  if (!is.null(control$method)) if (control$method != "quadra") ctrl$threshold <- 1e-2
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$method <- switch(ctrl$method, quadra = "QUADRA", pathwise = "PATHWISE", fista = "FISTA", 0)
  ctrl$normalize <- normalize
  ctrl$beta0  <- beta0
  
  ## ============================================
  ## INSTANTIATE THE DATA MODEL
  ##
  myData <- DataModel$new(
    covariates  = x,
    outcome     = y,
    cov_struct  = struct,
    obs_weights = weights
  )
  
  ## ============================================
  ## INSTANTIATE THE PENALIZED MODEL
  ##
  myModel <- LavaFit$new(
    data      = myData,
    intercept = intercept,
    regParam  = list(lambda = lambda1, 
                     gamma  = lambda2,
                     eta = 1,
                     lambda_factor = penscale, 
                     min_ratio = minratio, n_lambda = nlambda1)
  )
  
  ## ============================================
  ## FIT THE MODEL WITH ACTIVE SET ALGORITHM
  ##
  if (ctrl$verbose > 0) cat("\nModel fitting and optimization")
  myModel$fit(ctrl)

  ## ============================================
  ## POSTREATMENT + SEND BACK THE RESULTING MODEL
  ##
  if (ctrl$verbose > 0) cat("\nPost-treatment")
  myModel$debias <- refit
  myModel$criteria()
  myModel
}
