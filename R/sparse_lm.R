#' Fit a linear model with sparse regularization
#'
#' Adjust a linear model with sparse regularization.
#' We also add a (possibly structured) \eqn{\ell_2}{l2}-norm
#' (ridge-like). The solution path is computed at a grid of values for the
#' \eqn{\ell_1}{l1}-penalty, fixing the amount of \eqn{\ell_2}{l2}
#' regularization. See details for the criterion optimized.
#'
#' @param x matrix of features, possibly sparsely encoded
#' (experimental). Do NOT include intercept. When normalized is
#' `TRUE`, coefficients will then be rescaled to the original
#' scale.
#'
#' @param y response vector.
#'
#' @param lambda1 sequence of decreasing \eqn{\ell_1}{l1}-penalty
#' levels. If `NULL` (the default), a vector is generated with
#' `nlambda1` entries, starting from a guessed level
#' `lambda1.max` where only the intercept is included, then
#' shrunken to `minratio*lambda1.max`.
#'
#' @param lambda2 real scalar; tunes the \eqn{\ell_2}{l2} penalty in
#' the Elastic-net. Default is 0.01. Set to 0 to recover the Lasso.
#'
#' @param weights vector with real positive values that weight the
#' observations (like in weighted least square).
#' Default sets all weights to 1.
#'
#' @param penscale vector with real positive values that weight the
#' penalty of each feature. Default sets all weights to 1.
#'
#' @param struct matrix structuring the coefficients, possibly
#' sparsely encoded. Must be at least positive semidefinite (this is
#' checked internally). If `NULL` (the default), the identity matrix is
#' used. See details below.
#'
#' @param intercept logical; indicates if an intercept should be
#' included in the model. Default is `TRUE`.
#'
#' @param normalize logical; indicates if variables should be
#' normalized to have unit L2 norm before fitting.  Default is
#' `TRUE`.
#'
#' @param refit logical: indicates if the non null coefficients should be
#' refit to avoid excessive  bias. Default is FALSE. Can be changed later
#' (both raw and refit coefficients are stored).
#'
#' @param nlambda1 integer that indicates the number of values to put
#' in the `lambda1` vector.  Ignored if `lambda1` is
#' provided.
#'
#' @param minratio minimal value of \eqn{\ell_1}{l1}-part of the
#' penalty that will be tried, as a fraction of the maximal
#' `lambda1` value. A too small value might lead to instability
#' at the end of the solution path corresponding to small
#' `lambda1` combined with \eqn{\lambda_2=0}{lambda2=0}.  The
#' default value tries to avoid this, adapting to the
#' '\eqn{n<p}{n<p}' context. Ignored if `lambda1` is provided.
#'
#' @param maxfeat integer; limits the number of features ever to
#' enter the model; i.e., non-zero coefficients for the Elastic-net:
#' the algorithm stops if this number is exceeded and `lambda1`
#' is cut at the corresponding level. Default is
#' `min(nrow(x),ncol(x))` for small `lambda2` (<0.01) and
#' `min(4*nrow(x),ncol(x))` otherwise. Use with care, as it
#' considerably changes the computation time.
#'
#' @param beta0 a starting point for the vector of parameter. By default,
#' will initialized zero. May save time in some situation.
#'
#' @param control list of argument controlling low level options of
#' the algorithm --use with care and at your own risk-- :
#' * `verbose`: integer; activate verbose mode --this one is not
#' too risky!-- set to `0` for no output; `1` for warnings only,
#' and `2` for tracing the whole progression. Default is `1`.
#' Automatically set to `0` when the method is embedded within
#' cross-validation or stability selection.
#' * `timer`: logical; use to record the timing of the
#' algorithm. Default is `FALSE`.
#' * `maxiter` the maximal number of iteration used in the active set algorithm
#' to solve the problem for a given value of lambda1 . Default is 50.
#' * `method` a string for the underlying solver used. Either
#' `"quadra"`, `"fista"` or `"pgd"`. Default is `"quadra"`.
#' * `factmat` Boolean indicating if matrix factorization should be used to
#' solve the sub-system. If `TRUE` (the default), a Cholesky decomposition is
#' maintained along the path. If `FALSE`, the sub-system are solved with
#' a conjugate gradient algorithm.
#' * `threshold` a threshold for convergence. The algorithm stops
#' when the optimality conditions are fulfill up to this threshold.
#' Default is `1e-6`.
#' * `monitor` indicates if a monitoring of the convergence should be
#' recorded, by computing a lower bound between the current solution and
#' the optimum: when `'0'` (the default), no monitoring is provided;
#' when `'1'`, the bound derived in Grandvalet et al. is computed; when
#' `'>1'`, the Fenchel duality gap is computed along the algorithm.
#'
#' @param type string indicating the sparse variant to be fitted.
#' Could be "l1", "mcp" or "scad". Default is "l1".
#' be careful as scad and mcp are still experimental and have not been fully tested yet
#'
#' @param eta real positive scalar for tuning SCAD or MCP penalties.
#' Default is 3.7. Ignored when type == "l1".
#'
#' @details The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda_1,\lambda_2} = \arg \min_{\beta} \frac{1}{2}
#' (y - X \beta)^T (y - X \beta) + \lambda_1 pen_{\eta}(D \beta) +
#' \frac{\lambda_2}{2} \beta^T S \beta, }} \if{html}{\out{ 
#' &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
#' argmin<sub>&beta;</sub> 1/2 RSS(&beta;) + &lambda;<sub>1</sub>
#' pen<sub>&eta;</sub>(D &beta;) + &lambda;/2 <sub>2</sub>
#' &beta;<sup>T</sup> S &beta;,  }}
#' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
#' RSS(beta) + lambda1 |D beta|1 + lambda2 beta' S beta,}} where
#' \eqn{D}{D} is a diagonal matrix, whose diagonal terms are provided
#' as a vector by the \code{penscale} argument. The \eqn{\ell_2}{l2}
#' structuring matrix \eqn{S}{S} is provided via the `struct`
#' argument, a positive semidefinite matrix (possibly of class
#' `Matrix`).
#'
#' @return an object with class [SparseFit], inheriting from [QuadrupenFit].
#'
#' @seealso See also [SparseFit]
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
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- "relevant"
#'
#' ## Lasso
#' plot(lasso(x, y), label=labels)
#'
#' ## SCAD
#' plot(scad(x, y), label=labels)
#'
#' ## MCP
#' plot(mcp(x, y), label=labels)
#'
#' ## Elastic-net
#' plot(elastic_net(x,y,lambda2=1), label=labels)
#'
#' ## Structured Elastic-net (l2-structuring prior)
#' plot(elastic_net(x,y,lambda2=3,struct=solve(Sigma)), label=labels)
#'
#' ## SCAD + L2
#' plot(scad(x,y, eta = 3.7, lambda2=1), label=labels)
#'
#' ## MCP + L2
#' plot(mcp(x, y, eta = 3, lambda2=1), label=labels)
#'
#' @export
sparse_lm <- function(x,
                      y,
                      type      = c("l1", "mcp", "scad"),
                      lambda1   = NULL,
                      lambda2   = 0.01,
                      eta       = 3.7,
                      weights   = rep(1,nrow(x)),
                      penscale  = rep(1,ncol(x)),
                      struct    = Matrix::Diagonal(ncol(x), 1),
                      intercept = TRUE,
                      normalize = TRUE,
                      refit     = FALSE,
                      nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                      minratio  = ifelse(nrow(x) <= ncol(x), 1e-2, 1e-4),
                      maxfeat   = ifelse(lambda2 < 1e-2, min(nrow(x),ncol(x)), min(4*nrow(x),ncol(x))),
                      beta0     = numeric(ncol(x)),
                      control   = list()) {

  type <- match.arg(type)
  if (type == "mcp") stopifnot(eta > 1)
  if (type == "scad") stopifnot(eta > 2)

  ## ============================================
  ## RECOVER LOW LEVEL CONFIGURATION
  ##
  ctrl <- optim_enet_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$method <- switch(ctrl$method, quadra = "QUADRA", fista = "FISTA", pgd = "PGD", 0)
  ctrl$factmat <- ctrl$method == "QUADRA"
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
  myModel <- SparseFit$new(
    data      = myData,
    intercept = intercept,
    type      = type,
    regParam  = list(lambda = lambda1,
                     gamma = lambda2,
                     eta = eta,
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
