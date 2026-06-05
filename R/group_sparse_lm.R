#' Fit a linear model with (sparse) group regularisation
#'
#' Adjust a linear model with (sparse) group regularization, that is, a
#' mixture of an element-wise \eqn{\ell_1}{l1} norm and a group-wise mixed-norm
#' (either \eqn{\ell_1/\ell_2}{l1/l2}, \eqn{\ell_1/\ell_\infty}{l1/linf} or
#' cooperative). We also add a (possibly structured) \eqn{\ell_2}{l2}-norm (ridge-like).
#' The solution path is computed on an automatically tuned  grid of values for
#' the sparse group penalty. The mixture coefficient and the amount of ridge-like
#' regularization are fixed by the user.  See details for the criterion optimized.
#'
#' @inheritParams sparse_lm
#'
#' @param alpha real scalar in (0,1); tunes mixture between \eqn{\ell_1}{l1}
#' group penalties. Default is 0.0 (standard group-lasso).
#'
#' @param group vector of integers indicating group belonging. Must
#' match the number of column in \code{x}. Must be SORTED integers
#' starting from 1.
#'
#' @param type string indicating the sparse-group variant to be fitted.
#' Could be "l2", "coop", or "linf". Default is "l2" (regular Group-Lasso)
#'
#' @details The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda_1,\lambda_2} = \arg \min_{\beta} \frac{1}{2}
#' (y - X \beta)^T (y - X \beta) + \lambda_1 \left[ (1-\alpha)\,
#' \Omega_g(\beta) + \alpha \| D \beta \|_1 \right] +
#' \frac{\lambda_2}{2} \beta^T S \beta, }}
#' \if{html}{\out{ 
#' &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
#' argmin<sub>&beta;</sub> 1/2 RSS(&beta;) + &lambda;<sub>1</sub>
#' [(1-&alpha;) &Omega;<sub>g</sub>(&beta;) + &alpha; &#124; D &beta; &#124;<sub>1</sub>]
#' + &lambda;<sub>2</sub>/2 &beta;<sup>T</sup> S &beta;,
#'  }}
#' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
#' RSS(beta) + lambda1 [(1-alpha) Omega_g(beta) + alpha |D beta|_1]
#' + lambda2/2 beta' S beta,}}
#' where \eqn{\Omega_g}{Omega_g} is the group-wise mixed norm:
#' \eqn{\ell_1/\ell_2}{l1/l2} (Group-Lasso),
#' \eqn{\ell_1/\ell_\infty}{l1/linf}, or cooperative (Clime), controlled
#' by the \code{type} argument; \eqn{D}{D} is a diagonal matrix whose
#' diagonal terms are given by \code{penscale}; \eqn{\alpha}{alpha}
#' tunes the mixture between the group and element-wise penalties;
#' and \eqn{S}{S} is the \eqn{\ell_2}{l2} structuring matrix provided
#' via \code{struct}, a positive semidefinite matrix (possibly of
#' class \code{Matrix}).
#'
#' @return an object with class [SparseGroupFit], inheriting from [QuadrupenFit].
#'
#' @seealso See also [QuadrupenFit]
#'
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
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
#' ## Various sparse group linear models without/with an additional l2 regularization term
#' ## and with structuring prior
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- "relevant"
#'
#' ## Group-Lasso
#' plot(group_lasso(x, y, grp), label=labels)
#'
#' ## Sparse Group-Lasso
#' plot(sparse_group_lasso(x, y, grp, alpha = 0.75), label=labels)
#'
#' \dontrun{
#' 
#' ## Sparse Group-Lasso + L2 regularization
#' plot(group_sparse_lm(x, y, grp, type = "l2", alpha = .75, lambda2=.5),
#'  label=labels)
#' plot(group_sparse_lm(x, y, grp, type = "l2", alpha = .75, lambda2=10),
#'  label=labels)
#' plot(group_sparse_lm(x, y, grp, type = "l2", alpha = .75, lambda2=10,
#'  struct=solve(Sigma)), label=labels)
#'
#' ## Group-Lasso L1/LINF
#' plot(group_l1linf(x, y, grp), label=labels)
#'
#' ## Sparse Group-Lasso L1/LINF
#' plot(sparse_group_l1linf(x, y, grp, alpha = 0.75), label=labels)
#'
#' ## Sparse L1/LINF Group-Lasso + L2 regularization
#' plot(group_sparse_lm(x, y, grp, type = "linf", alpha = .75, lambda2=.5),
#'   label=labels)
#' plot(group_sparse_lm(x, y, grp, type = "linf", alpha = .75, lambda2=10),
#'   label=labels)
#' plot(group_sparse_lm(x, y, grp, type = "linf", alpha = .75, lambda2=10,
#'   struct=solve(Sigma)), label=labels)
#'
#' ## Cooperative-Lasso
#' plot(coop_lasso(x, y, grp), label=labels)
#'
#' ## Sparse Cooperative-Lasso
#' plot(sparse_coop_lasso(x, y, grp, alpha = 0.75), label=labels)
#'
#' ## Sparse Cooperative-Lasso + L2 regularization
#' plot(group_sparse_lm(x, y, grp, type = "coop", alpha = .75, lambda2=.5),
#'  label=labels)
#' plot(group_sparse_lm(x, y, grp, type = "coop", alpha = .75, lambda2=10),
#'  label=labels)
#' plot(group_sparse_lm(x, y, grp, type = "coop", alpha = .75, lambda2=10,
#'  struct=solve(Sigma)), label=labels)
#' }
#' @export
group_sparse_lm <-
  function(x,
           y,
           group,
           type      = c("l2", "coop", "linf"),
           lambda1   = NULL,
           lambda2   = 0.01,
           alpha     = 0.0,
           weights   = rep(1,nrow(x)),
           penscale  = sqrt(tabulate(group)),
           struct    = Matrix::Diagonal(ncol(x), 1),
           intercept = TRUE,
           normalize = TRUE,
           refit     = FALSE,
           nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
           minratio  = 1e-2,
           maxfeat   = ifelse(lambda2 < 1e-2, min(2*nrow(x),ncol(x)), min(4*nrow(x),ncol(x))),
           beta0     = numeric(ncol(x)),
           control   = list()) {

    stopifnot(alpha < 1 && alpha >= 0)
    stopifnot(!is.unsorted(group))
    type <- match.arg(type)

    ## ============================================
    ## RECOVER LOW LEVEL CONFIGURATION
    ##
    ctrl <- optim_grp_default(ncol(x))
    if (is.null(control$method) && type == "l2" && alpha == 0.0) control$method <- "quadra"
    ctrl[names(control)] <- control # default overwritten by user specifications
    ctrl$method  <- switch(ctrl$method, quadra = "QUADRA", fista = "FISTA", pgd = "PGD", 0)
    ctrl$factmat <- ctrl$method == "QUADRA" && type == "l2"
    ctrl$normalize <- normalize
    ctrl$beta0 <- beta0
    ctrl$maxfeat <- maxfeat

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
    myModel <- SparseGroupFit$new(
      data      = myData,
      intercept = intercept,
      group     = group,
      type      = type,
      regParam  = list(lambda = lambda1,
                       gamma  = lambda2,
                       alpha  = alpha,
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
