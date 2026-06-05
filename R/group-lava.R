#' Fit a linear model with group-lava regularization
#'
#' Adjust a the group-lava regularized linear models, that is a lava transformation
#' of the data plus a mixture of either a (possibly weighted)
#' \eqn{\ell_1/\ell_2}{l1/l2}- or
#' \eqn{\ell_1/\ell_\infty}{l1/linf}-norm, and a (possibly
#' structured) \eqn{\ell_2}{l2}-norm (ridge-like). The solution path
#' is computed at a grid of values for the
#' \eqn{\ell_1/\ell_q}{l1/lq}-penalty. See details for the criterion
#' optimized.
#'
#' @inheritParams elastic_net
#'
#' @param group vector of integers indicating group belonging. Must
#' match the number of column in \code{x}. Must be SORTED integers
#' starting from 1.
#' 
#' @param type string indicating whether the \eqn{\ell_1/\ell_2}{l1/l2} or the
#' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted. Could be "linf" or 
#' "l2", default is "l2"
#'
#' @details The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\theta}_{\lambda_1,\lambda_2} = \arg \min_{\theta = \beta + \delta} \frac{1}{2 n} (y - X
#' (\beta + \delta))^T (y - X (\beta + \delta)) + \lambda_1 \Omega_g(\beta) + \frac{\lambda_2}{2} \delta^T S \delta, }}
#' \if{html}{\out{  &theta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
#' argmin<sub>&theta; = &beta;+&delta;</sub> 1/2n RSS(&beta; + &delta;) + &lambda;<sub>1</sub>
#' &Omega;<sub>g</sub>(&beta;) + &lambda;<sub>2</sub>/2 &delta;<sup>T</sup> S
#' &delta;,  }}
#' \if{text}{\deqn{theta.hat(lambda1,lambda2) = argmin_{theta = beta + delta} 1/2n
#' RSS(beta + delta) + lambda1 Omega_g(beta) + lambda2/2 delta' S delta,}}
#' where \eqn{\Omega_g}{Omega_g} is the group-wise mixed norm:
#' \eqn{\ell_1/\ell_2}{l1/l2} (Group-LAVA) or
#' \eqn{\ell_1/\ell_\infty}{l1/linf}, controlled by the \code{type} argument.
#' The \eqn{\ell_2}{l2} structuring matrix \eqn{S}{S} is provided via \code{struct}.
#'
#' @return an object with class [GroupLavaFit], inheriting from [QuadrupenFit].
#' 
#' @references Chernozhukov, Victor, Christian Hansen, and Yuan Liao. "A lava attack on the recovery of sums of 
#' dense and sparse signals." The Annals of Statistics (2017): 39-76. <doi:10.1214/16-AOS1434>
#' 
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta  <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' delta <- runif(sum(c(25,10,25,10,25)),-.1,.1)
#' grp  <- rep(1:5, c(25,10,25,10,25)) 
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
#' 
#' \dontrun{
#' ## Standard Group-Lasso
#' plot(group_lava(x,y,grp), label=labels)
#' plot(group_lava(x,y,grp, lambda2=.5), label=labels)
#' plot(group_lava(x,y,grp, lambda2=10), label=labels)
#' plot(group_lava(x,y,grp, lambda2=10,struct=solve(Sigma)), label=labels)
#' 
#' ## L1/LINF Group-Lasso
#' plot(group_lava(x, y, grp, type = "linf"), label=labels)
#' plot(group_lava(x, y, grp, type = "linf", lambda2=.5), label=labels)
#' plot(group_lava(x, y, grp, type = "linf", lambda2=10), label=labels)
#' plot(group_lava(x, y, grp, type = "linf", lambda2=10,struct=solve(Sigma)), label=labels)
#' 
#' ## Cooperative-Lasso
#' plot(group_lava(x, y, grp, type = "coop"), label=labels)
#' plot(group_lava(x, y, grp, type = "coop", lambda2=.5), label=labels)
#' plot(group_lava(x, y, grp, type = "coop", lambda2=10), label=labels)
#' plot(group_lava(x, y, grp, type = "coop", lambda2=10,struct=solve(Sigma)), label=labels)
#' }
#'
#' @export
group_lava <- function(x,
                       y,
                       group,
                       type        = c("l2", "coop", "linf"),
                       lambda1   = NULL,
                       lambda2   = 1,
                       weights   = rep(1, nrow(x)),
                       penscale  = sqrt(tabulate(group)),
                       struct    = Matrix::Diagonal(ncol(x), 1),
                       intercept = TRUE,
                       normalize = TRUE,
                       refit     = FALSE,
                       nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                       minratio  = 1e-2,
                       maxfeat   = ifelse(lambda2 < 1e-2, min(nrow(x),ncol(x)), min(4*nrow(x),ncol(x))),
                       beta0     = numeric(ncol(x)),
                       control   = list()) {

  ## ============================================
  ## RECOVER LOW LEVEL CONFIGURATION
  ##
  ctrl <- optim_grp_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  if (!is.null(control$method)) if (control$method != "quadra") ctrl$threshold <- 1e-2
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$method  <- switch(ctrl$method, quadra = "QUADRA", fista = "FISTA", 0)
  ctrl$factmat <- FALSE
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
  myModel <- GroupLavaFit$new(
    data      = myData,
    intercept = intercept,
    group     = group,
    type      = match.arg(type),
    regParam  = list(lambda = lambda1,
                     gamma  = lambda2,
                     alpha  = 0,
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
