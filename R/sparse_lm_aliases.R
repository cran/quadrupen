#' @rdname sparse_lm
#' @importFrom lifecycle badge deprecate_warn
#' @export
elastic.net <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.5,
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
                        control   = list(method = "quadra")) {

  lifecycle::deprecate_warn("1.1.0", "elastic.net()", "elastic_net()")

  out <- sparse_lm(x,
                   y,
                   type      = "l1",
                   lambda1   = lambda1,
                   lambda2   = lambda2,
                   eta       = 0,
                   weights   = weights,
                   penscale  = penscale,
                   struct    = struct,
                   intercept = intercept,
                   normalize = normalize,
                   refit     = refit,
                   nlambda1  = nlambda1,
                   minratio  = minratio,
                   maxfeat   = maxfeat,
                   beta0     = beta0,
                   control   = control)
  out
}

#' @rdname sparse_lm
#' @export
elastic_net <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.5,
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
                        control   = list(method = "quadra")) {

  out <- sparse_lm(x,
                   y,
                   type      = "l1",
                   lambda1   = lambda1,
                   lambda2   = lambda2,
                   eta       = 0,
                   weights   = weights,
                   penscale  = penscale,
                   struct    = struct,
                   intercept = intercept,
                   normalize = normalize,
                   refit     = refit,
                   nlambda1  = nlambda1,
                   minratio  = minratio,
                   maxfeat   = maxfeat,
                   beta0     = beta0,
                   control   = control)
  out
}

#' @rdname sparse_lm
#' @export
lasso <- function(x,
                  y,
                  lambda1   = NULL,
                  weights   = rep(1,nrow(x)),
                  penscale  = rep(1,ncol(x)),
                  intercept = TRUE,
                  normalize = TRUE,
                  refit     = FALSE,
                  nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                  minratio  = ifelse(nrow(x) <= ncol(x), 1e-2, 1e-4),
                  maxfeat   = min(nrow(x),ncol(x)),
                  beta0     = numeric(ncol(x)),
                  control   = list(method = "quadra")) {

  out <- sparse_lm(x,
                   y,
                   type = "l1",
                   lambda1   = lambda1,
                   lambda2   = 0,
                   eta       = 0,
                   weights   = weights,
                   penscale  = penscale,
                   intercept = intercept,
                   normalize = normalize,
                   refit     = refit,
                   nlambda1  = nlambda1,
                   minratio  = minratio,
                   maxfeat   = maxfeat,
                   beta0     = beta0,
                   control   = control)
  out
}

#' @rdname sparse_lm
#' @export
mcp <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.0,
                        eta       = 3,
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
                        control   = list(method = "quadra")) {

  out <- sparse_lm(x,
                   y,
                   type      = "mcp",
                   lambda1   = lambda1,
                   lambda2   = lambda2,
                   eta       = eta,
                   weights   = weights,
                   penscale  = penscale,
                   struct    = struct,
                   intercept = intercept,
                   normalize = normalize,
                   refit     = refit,
                   nlambda1  = nlambda1,
                   minratio  = minratio,
                   maxfeat   = maxfeat,
                   beta0     = beta0,
                   control   = control)
  out
}

#' @rdname sparse_lm
#' @export
scad <- function(x,
                 y,
                 lambda1   = NULL,
                 lambda2   = 0.0,
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
                 control   = list(method = "quadra")) {

  out <- sparse_lm(x,
                   y,
                   type      = "scad",
                   lambda1   = lambda1,
                   lambda2   = lambda2,
                   eta       = eta,
                   weights   = weights,
                   penscale  = penscale,
                   struct    = struct,
                   intercept = intercept,
                   normalize = normalize,
                   refit     = refit,
                   nlambda1  = nlambda1,
                   minratio  = minratio,
                   maxfeat   = maxfeat,
                   beta0     = beta0,
                   control   = control)
  out
}
