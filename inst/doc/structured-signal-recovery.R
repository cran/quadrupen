## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup, message = FALSE---------------------------------------------------
library(quadrupen)
library(Matrix)
library(ggplot2)

## ----simulation---------------------------------------------------------------
set.seed(42)
p    <- 95
n    <- 50
beta <- rep(c(0, 1, 0, -1, 0), c(25, 10, 25, 10, 25))

## Block-correlated predictors: high within-segment, zero across segments
rho <- 0.75
Soo <- toeplitz(rho^(0:24))      # Toeplitz correlation within zero segments
Sww <- matrix(rho, 10, 10)       # constant correlation within active segments
Sigma <- bdiag(Soo, Sww, Soo, Sww, Soo)
diag(Sigma) <- 1

X <- as.matrix(matrix(rnorm(p * n), n, p) %*% chol(Sigma))
y <- X %*% beta + rnorm(n, 0, 10)

## Segment labels (for plot legends)
segments      <- rep(1:5, c(25, 10, 25, 10, 25))
segment_names <- paste0("seg", 1:5)
seg_labels    <- segment_names[segments]

## ----plot-beta, fig.height = 3------------------------------------------------
data.frame(index = seq_along(beta), beta = beta, segment = seg_labels) |>
  ggplot(aes(x = index, y = beta, colour = segment)) +
  geom_step() + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  labs(x = "Predictor index", y = expression(beta), title = "True coefficient vector") +
  theme_bw()

## ----graph--------------------------------------------------------------------
## Adjacency matrix: clique within each segment, no edges across
Ioo <- matrix(1, 25, 25)
Iww <- matrix(1, 10, 10)
A   <- as.matrix(bdiag(Ioo, Iww, Ioo, Iww, Ioo))
diag(A) <- 0

## Graph Laplacian: L = D - A  (regularized for positive definiteness)
L2        <- -A
diag(L2)  <- colSums(A) + 1e-2

## ----graph-viz, fig.height = 4, fig.width = 5, fig.align = "center"-----------
if (requireNamespace("igraph", quietly = TRUE)) {
  g  <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
  vertex_colors <- cols[segments]
  igraph::plot.igraph(
    g, vertex.color = vertex_colors, vertex.size = 3,
    vertex.label = NA, edge.width = 0.1,
    main = "Community graph (colour = segment)"
  )
  legend("topleft", legend = segment_names, fill = cols, bty = "n", cex = 0.8)
}

## ----lasso--------------------------------------------------------------------
fit_lasso <- lasso(X, y, intercept = FALSE)
fit_lasso

## ----lasso-path---------------------------------------------------------------
fit_lasso$plot_path(labels = seg_labels)

## ----fused-chain--------------------------------------------------------------
fit_fused_chain <- fused_lasso(X, y, lambda2 = 5, intercept = FALSE)
fit_fused_chain

## ----fused-chain-path---------------------------------------------------------
fit_fused_chain$plot_path(labels = seg_labels)

## ----fused-community----------------------------------------------------------
fit_fused_comm <- fused_lasso(X, y, lambda2 = 5, struct = A, intercept = FALSE)
fit_fused_comm

## ----fused-community-path-----------------------------------------------------
fit_fused_comm$plot_path(labels = seg_labels)

## ----ridge-std----------------------------------------------------------------
fit_ridge_std <- ridge(X, y, intercept = FALSE)
fit_ridge_std$plot_path(labels = seg_labels)

## ----ridge-struct-------------------------------------------------------------
fit_ridge_struct <- ridge(X, y, struct = L2, lambda_max = 1000, intercept = FALSE)
fit_ridge_struct$plot_path(labels = seg_labels)

## ----enet-struct--------------------------------------------------------------
fit_enet_struct <- elastic_net(X, y, lambda2 = 5, struct = L2, intercept = FALSE)
fit_enet_struct

## ----enet-struct-path---------------------------------------------------------
fit_enet_struct$plot_path(labels = seg_labels)

## ----enet-debias--------------------------------------------------------------
fit_enet_struct$debias <- TRUE
fit_enet_struct$plot_path(labels = seg_labels)
fit_enet_struct$debias <- FALSE

## ----recovery-----------------------------------------------------------------
methods <- list(
  "Lasso"                  = fit_lasso,
  "Fused Lasso (chain)"    = fit_fused_chain,
  "Fused Lasso (community)"= fit_fused_comm,
  "Ridge (standard)"       = fit_ridge_std,
  "Ridge (structured)"     = fit_ridge_struct,
  "Elastic-net (structured)"= fit_enet_struct
)

extract_coef <- function(fit) {
  b <- fit$get_model("BIC")
  if ("Intercept" %in% names(b)) b <- b[-1]   # drop intercept
  b
}

df_coef <- do.call(rbind, lapply(names(methods), function(nm) {
  b <- extract_coef(methods[[nm]])
  data.frame(
    method  = factor(nm, levels = names(methods)),
    index   = seq_along(b),
    value   = as.numeric(b),
    segment = seg_labels
  )
}))

df_true <- data.frame(
  index   = seq_along(beta),
  beta    = beta,
  segment = seg_labels
)

## ----recovery-plot, fig.height = 9, fig.width = 7-----------------------------
## Background shading for each segment
seg_bounds <- data.frame(
  xmin = c(1, 26, 36, 61, 71),
  xmax = c(25, 35, 60, 70, 95),
  fill = factor(1:5)
)

ggplot(df_coef, aes(x = index, y = value)) +
  geom_rect(data = seg_bounds,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.08, inherit.aes = FALSE) +
  scale_fill_manual(values = c("#E69F00","#56B4E9","#009E73","#F0E442","#CC79A7"),
                    guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_step(data = df_true, aes(y = beta), colour = "black", linewidth = 0.8,
            linetype = "dashed") +
  geom_point(aes(colour = segment), size = 0.8, alpha = 0.7) +
  scale_colour_manual(values = c("#E69F00","#56B4E9","#009E73","#F0E442","#CC79A7")) +
  facet_wrap(~ method, ncol = 2) +
  labs(x = "Predictor index", y = "Estimated coefficient",
       title = "Coefficient recovery by method",
       subtitle = "Dashed line: true β",
       colour = "Segment") +
  theme_bw() +
  theme(legend.position = "bottom")

