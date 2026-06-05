#' A function for fitting generalized fused-Lasso problems
#' 
#' This function fits the standard version of the fused lasso. 
#' It can take a general matrix x and allows for possible weights on the 
#' `lambda1` and `lambda2` penalties. 
#' 
#' @inheritParams elastic_net
#' 
#' @param struct Description of the graph that corresponds to the `lambda2` penalty structure. 
#' If \code{NULL} (the default) a chain graph is assumed, like in the standard fused-lasso.
#' If a matrix is given, interpreted as 
#' a symmetric adjacency matrix
#' @param pen_fused penalty used for fusing the variables (either L1, L2 or Huber). Default is L1 
#' 
#' @param control list of argument controlling low level options of
#' the algorithm:
#' * `verbose`: logical; verbose mode
#' * `timer`: logical; use to record the timing of the
#' algorithm. Default is `FALSE`.
#' * `maxiterin` Maximum number of iterations in the inner loop to run.
#' * `maxiterout` Maximum number of iterations in the outer loop to run.
#' * `maxactivation` Maximum number of previously inactive variables to activate at the same time
#' * `accuracy` Accuracy at which the algorithm will stop.
#' * `fusioncheck` Should the fused sets be checked for separation?
#' * `verbose` Should the function give some output what it is doing?
#' 
#' @return an object with class [FusedLassoFit], inheriting from [QuadrupenFit].
#' 
#' @details The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda_1,\lambda_2} = \arg \min_{\beta} \frac{1}{2}
#' (y - X \beta)^T (y - X \beta) + \lambda_1 \|w \circ \beta \|_{1} +
#' \frac{\lambda_2}{2} \sum_{i \sim j} w_{ij} |\beta_i - \beta_j|, }} \if{html}{\out{ 
#' &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
#' argmin<sub>&beta;</sub> 1/2 RSS(&beta;) + &lambda;<sub>1</sub>
#' &#124; D &beta; &#124;<sub>1</sub> + &lambda;/2 <sub>2</sub>
#' &beta;<sup>T</sup> S &beta;,  }}
#' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
#' RSS(beta) + lambda1 |D beta|1 + lambda2 sum{(i,j) in G} w(ij) |beta_j -beta_i|,}} where
#' \eqn{D}{D} is a diagonal matrix, whose diagonal terms are provided
#' as a vector by the \code{penscale} argument. The \eqn{\ell_1}{l1} fusion penalty 
#' is structured by a possibly weighted graph \eqn{G}{G} provided via the `struct`
#' argument, as a symmetric (undirected) adjacency matrix.
#' 
#' @author Original code by Holger Hoefling, refactoring by Julien Chiquet
#' 
#' @examples 
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
#' res <- fused_lasso(x, y, lambda2=5)
#' G <- igraph::make_ring(ncol(x)) |> igraph::as_adjacency_matrix(sparse = FALSE)
#' resG <- fused_lasso(x, y, lambda2=5, struct = G)
#' plot(res)
#' plot(resG)
#' 
#' @importFrom methods as
#' @export 
fused_lasso <- function(
    x, 
    y, 
    lambda1   = NULL,                       
    lambda2   = 1,
    pen_fused  =  c("L1", "L2", "Huber"),
    penscale  = rep(1,ncol(x)),
    struct    = NULL,
    intercept = TRUE,
    normalize = TRUE,
    nlambda1  = ifelse(is.null(lambda1),50,length(lambda1)),
    minratio  = 1e-2,
    maxfeat   = ifelse(lambda2 < 1, min(nrow(x),ncol(x)), min(2*nrow(x),ncol(x))),
    beta0     = rep(0, ncol(x)),
    control   = list()) {
  
  ## ============================================
  ## RECOVER LOW LEVEL CONFIGURATION
  ##
  ctrl <- optim_fused_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$beta0  <- beta0
  ctrl$mu0    <- 0.001
  ctrl$pen_fused  <- match.arg(pen_fused)
  ctrl$normalize <- normalize
  ctrl$fusioncheck <- 
    switch(ctrl$fusioncheck, "all" = 2L, "active" = 1L, "none" = 0L, "naive" = -1L, -2L)
  if (ctrl$pen_fused == "L2")  ctrl$fusioncheck <- 0
  if (is.null(struct))   struct <- chain_graph(ncol(x))
  if (is.matrix(struct)) struct <- conn_from_adj(struct)
  
  ## ============================================
  ## INSTANTIATE THE DATA MODEL
  ##
  if (ctrl$verbose > 0) cat ("\nData pretreatment")
  myData <- DataModel$new(
    covariates  = as(x, "dgCMatrix"),
    outcome     = y,
    cov_struct  = struct
  )

  ## ============================================
  ## INSTANTIATE THE PENALIZED MODEL
  ##
  myModel <- FusedLassoFit$new(
    data      = myData,
    intercept = intercept,
    regParam  = list(l1 = lambda1, l2 = lambda2, 
                     l1_weights = penscale, 
                     min_ratio = minratio, n_lambda1 = nlambda1)
    
  )
  
  ## ============================================
  ## FIT THE MODEL WITH ACTIVE SET ALGORITHM
  ##
  if (ctrl$verbose > 0) cat("\nModel fitting and optimization")
  myModel$fit(ctrl)

  ## ===========================
  ## Post-Treatments
  ## 
  if (ctrl$verbose > 0) cat("\nPost-treatment")
  myModel$criteria()
  myModel
}
