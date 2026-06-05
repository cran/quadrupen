## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(quadrupen)
data("Birthwt", package = "grpreg")
y     <- Birthwt$bwt[-130]    ## outlier
X     <- Birthwt$X[-130, ]
group <- as.integer(Birthwt$group)[-130]

## ----data-prep----------------------------------------------------------------
group_names <- levels(Birthwt$group)
var_labels  <- group_names[group]
cat("Groups (", length(group_names), "):", paste(group_names, collapse = ", "), "\n")
cat("Group sizes:", tabulate(group), "\n")

## ----gl-fit-------------------------------------------------------------------
fit_gl <- group_lasso(X, y, group)
fit_gl

## ----cl-fit-------------------------------------------------------------------
fit_cl <- coop_lasso(X, y, group)
fit_cl

## ----sgl-fit------------------------------------------------------------------
fit_sgl <- sparse_group_lasso(X, y, group, alpha = 0.5)
fit_sgl

## ----path-plots, fig.show='hold', out.width='50%'-----------------------------
fit_gl$plot_path(xvar = "fraction", log_scale = FALSE, labels = var_labels)
fit_cl$plot_path(labels = var_labels)
fit_sgl$plot_path(labels = var_labels)
fit_sgl$plot_path(standardize = FALSE)

## ----criteria-----------------------------------------------------------------
fit_gl$criteria()

## ----criteria-plots, fig.show='hold', out.width='50%'-------------------------
fit_gl$plot(type = "criteria")

## ----criteria-plots-2, fig.show='hold', out.width='50%'-----------------------
fit_gl$information_criteria$plot(c("AIC", "BIC", "mBIC"))
fit_sgl$information_criteria$plot("GCV")

## ----get-model----------------------------------------------------------------
coef_gl  <- fit_gl$get_model("BIC")
coef_sgl <- fit_sgl$get_model("BIC")

active_gl  <- unique(group[coef_gl[-1]  != 0])
active_sgl <- unique(group[coef_sgl[-1] != 0])
cat("Group Lasso   — active groups (BIC):", group_names[active_gl],  "\n")
cat("Sparse GL     — active groups (BIC):", group_names[active_sgl], "\n")

## ----sgl-within---------------------------------------------------------------
active_vars_sgl <- which(coef_sgl[-1] != 0)
cat("Sparse GL — active predictors (BIC):", colnames(X)[active_vars_sgl], "\n")

## ----get-model-lambda---------------------------------------------------------
lambda_gl <- fit_gl$major_tuning
coef_gl_lambda <- fit_gl$get_model(lambda_gl[20])
cat("Non-zero Group Lasso coefficients for lambda =", round(lambda_gl[20], 3), "\n")
print(coef_gl_lambda[coef_gl_lambda != 0])

## ----cv, message=FALSE--------------------------------------------------------
set.seed(42)
fit_gl$cross_validate(K = 10, verbose = FALSE)
fit_sgl$cross_validate(K = 10, verbose = FALSE)

## ----cv-plots, fig.show='hold', out.width='50%'-------------------------------
fit_gl$plot(type = "crossval")
fit_sgl$plot(type = "crossval")

## ----cv-model-----------------------------------------------------------------
coef_gl_cv <- fit_gl$get_model("CV_min")
cat("Group Lasso — active groups (CV_min):",
    group_names[unique(group[coef_gl_cv[-1] != 0])], "\n")

## ----cv2D, message=FALSE------------------------------------------------------
set.seed(42)
fit_gl$cross_validate(K = 5, lambda2 = 10^seq(1, -3, len = 20), verbose = FALSE)

## ----cv-plots-2D, fig.show='hold', out.width='50%'----------------------------
fit_gl$plot(type = "crossval")
fit_gl$cross_validation$plotCV_1D(se = FALSE) + ggplot2::ylim(c(0.04, 0.08))

## ----predict------------------------------------------------------------------
y_hat_gl  <- fit_gl$predict(selection = "CV_min")
y_hat_cl  <- fit_cl$predict(selection = "BIC")
y_hat_sgl <- fit_sgl$predict(selection = "CV_min")

r2_gl  <- fit_gl$r_squared[fit_gl$get_model("CV_min", type = "index")]
r2_cl  <- fit_cl$r_squared[fit_cl$get_model("BIC",    type = "index")]
r2_sgl <- fit_sgl$r_squared[fit_sgl$get_model("CV_min", type = "index")]
cat(sprintf(
  "R² Group Lasso (CV_min): %.3f\nR² Coop Lasso  (BIC)   : %.3f\nR² Sparse GL   (CV_min): %.3f\n",
  r2_gl, r2_cl, r2_sgl
))

## ----stability, message=FALSE-------------------------------------------------
set.seed(42)
fit_gl$stability(n_subsamples = 200, verbose = FALSE)

## ----stability-plot-----------------------------------------------------------
fit_gl$plot(type = "stability", labels = var_labels)
fit_gl$stability_path$plot(nvarsel = 5, labels = var_labels)
colnames(X)[fit_gl$stability_path$selection(nvarsel = 5)]

## ----debias-------------------------------------------------------------------
fit_gl$debias <- TRUE
fit_gl$plot_path(labels = var_labels)

## ----debias-cv----------------------------------------------------------------
fit_gl$cross_validate(verbose = FALSE)
fit_gl$plot(type = "crossval")
fit_gl$debias <- FALSE
fit_gl$plot_path(labels = var_labels)

## ----no-debias----------------------------------------------------------------
fit_gl$debias <- FALSE
fit_gl$plot_path(labels = var_labels)

