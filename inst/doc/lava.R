## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup, message = FALSE---------------------------------------------------
library(quadrupen)
library(ggplot2)

## ----simulation---------------------------------------------------------------
set.seed(42)
n <- 100
p <- 50

## Sparse component: 5 large effects, the rest zero
beta  <- numeric(p)
beta[c(5, 15, 25, 35, 45)] <- c(3, -3, 3, -3, 3)

## Dense component: small effects on all predictors
delta <- rnorm(p, mean = 0, sd = 0.2)

## Total signal
theta <- beta + delta

## Predictors with moderate AR(1) correlation
rho   <- 0.5
Sigma <- toeplitz(rho^(0:(p - 1)))
X     <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
y     <- X %*% theta + rnorm(n, sd = 2)

## ----plot-truth, fig.height = 5-----------------------------------------------
df_truth <- data.frame(
  index     = rep(1:p, 3),
  value     = c(beta, delta, theta),
  component = factor(
    rep(c("Sparse (β)", "Dense (δ)", "Total (θ = β + δ)"), each = p),
    levels = c("Sparse (β)", "Dense (δ)", "Total (θ = β + δ)")
  )
)

ggplot(df_truth, aes(x = index, y = value, fill = component)) +
  geom_col() +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73"), guide = "none") +
  labs(x = "Predictor index", y = "Value",
       title = "True signal decomposition") +
  theme_bw()

## ----lava-fit-----------------------------------------------------------------
fit_lava <- lava(X, y, lambda2 = 1)
fit_lava

## ----paths, fig.show='hold', out.width='33%'----------------------------------
fit_lava$plot_path(component = "both",   labels = NULL)
fit_lava$plot_path(component = "sparse", labels = NULL)
fit_lava$plot_path(component = "dense",  labels = NULL)

## ----lava-criteria------------------------------------------------------------
fit_lava$information_criteria$plot(c("AIC", "BIC", "eBIC", "mBIC"))

## ----extraction---------------------------------------------------------------
idx <- fit_lava$get_model("BIC", type = "index")

theta_hat <- fit_lava$get_model("BIC")[-1]          # drop intercept
beta_hat  <- as.numeric(fit_lava$sparse_coef[, idx])
delta_hat <- as.numeric(fit_lava$dense_coef[, idx])

cat("Non-zero entries in sparse component (BIC):",
    sum(beta_hat != 0), "\n")
cat("Non-zero entries in dense component (BIC): ",
    sum(delta_hat != 0), "\n")
cat("Positions with |beta_hat| > 1:",
    which(abs(beta_hat) > 1), "\n")

## ----comparison-fit-----------------------------------------------------------
fit_lasso <- lasso(X, y)
fit_enet  <- elastic_net(X, y, lambda2 = 1)

## ----comparison-plot, fig.height = 6------------------------------------------
get_coef <- function(fit) {
  b <- fit$get_model("BIC")
  if ("Intercept" %in% names(b)) b <- b[-1]
  as.numeric(b)
}

method_levels <- c("True signal", "Lasso (BIC)", "Elastic-net (BIC)", "LAVA (BIC)")

## Non-LAVA panels: one bar per predictor
df_base <- rbind(
  data.frame(method = "True signal",       index = 1:p,
             value = theta,                fill_grp = "True signal"),
  data.frame(method = "Lasso (BIC)",       index = 1:p,
             value = get_coef(fit_lasso),  fill_grp = "Lasso"),
  data.frame(method = "Elastic-net (BIC)", index = 1:p,
             value = get_coef(fit_enet),   fill_grp = "Elastic-net")
)
df_base$method <- factor(df_base$method, levels = method_levels)

## LAVA panel: dense total drawn first, sparse component overlaid on top
df_lava_dense <- data.frame(
  method   = factor("LAVA (BIC)", levels = method_levels),
  index    = 1:p,
  value    = as.numeric(theta_hat),
  fill_grp = "LAVA dense (δ̂)"
)
df_lava_sparse <- data.frame(
  method   = factor("LAVA (BIC)", levels = method_levels),
  index    = 1:p,
  value    = beta_hat,
  fill_grp = "LAVA sparse (β̂)"
)

fill_palette <- c(
  "True signal"       = "grey20",
  "Lasso"             = "#56B4E9",
  "Elastic-net"       = "#E69F00",
  "LAVA dense (δ̂)"   = "#CC79A7",
  "LAVA sparse (β̂)"  = "#009E73"
)

ggplot(mapping = aes(x = index, y = value, fill = fill_grp)) +
  geom_col(data = df_base) +
  geom_col(data = df_lava_dense) +
  geom_col(data = df_lava_sparse) +
  facet_wrap(~ method, ncol = 2) +
  scale_fill_manual(
    values = fill_palette,
    breaks = c("LAVA sparse (β̂)", "LAVA dense (δ̂)"),
    name   = "LAVA decomposition"
  ) +
  labs(x = "Predictor index", y = "Estimated coefficient",
       title = "Coefficient recovery: Lasso vs. Elastic-net vs. LAVA") +
  theme_bw() +
  theme(legend.position = "bottom")

## ----cv, message = FALSE------------------------------------------------------
set.seed(42)
fit_lava$cross_validate(K = 10, verbose = FALSE)
fit_lava$plot(type = "crossval")

## ----cv-model-----------------------------------------------------------------
set.seed(42)
fit_lasso$cross_validate(K = 10, verbose = FALSE)
fit_enet$cross_validate(K = 10, verbose = FALSE)

idx_cv   <- fit_lava$get_model("CV_min", type = "index")
beta_cv  <- as.numeric(fit_lava$sparse_coef[, idx_cv])
cat("Non-zero entries in sparse component (CV_min):",
    sum(beta_cv != 0), "\n")
cat("Positions with |beta_cv| > 1:",
    which(abs(beta_cv) > 1), "\n")

## ----group-simulation---------------------------------------------------------
## 5 groups of 10 predictors; groups 2 and 4 are active in the sparse part
group  <- rep(1:5, each = 10)
beta_g <- rep(c(0, 2, 0, -2, 0), each = 10)
delta_g <- rnorm(p, mean = 0, sd = 0.2)
theta_g <- beta_g + delta_g

y2 <- X %*% theta_g + rnorm(n, sd = 2)

group_names <- paste0("grp", 1:5)
var_labels  <- group_names[group]

## ----group-lava-fit-----------------------------------------------------------
fit_gl <- group_lava(X, y2, group = group, lambda2 = 1)
fit_gl

## ----group-paths, fig.show='hold', out.width='50%'----------------------------
fit_gl$plot_path(component = "sparse", labels = var_labels)
fit_gl$plot_path(component = "dense",  labels = var_labels)

## ----group-extraction---------------------------------------------------------
idx_ebic   <- fit_gl$get_model("eBIC", type = "index")
beta_g_hat <- as.numeric(fit_gl$sparse_coef[, idx_ebic])

active_groups <- unique(group[beta_g_hat != 0])
cat("Active groups in sparse component (eBIC):", group_names[active_groups], "\n")

