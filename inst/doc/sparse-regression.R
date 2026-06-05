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
y <- Birthwt$bwt[-130] ## outlier
X <- Birthwt$X[-130, ]

## ----lasso-fit----------------------------------------------------------------
fit_lasso <- lasso(X, y, minratio = 1e-2)
fit_lasso

## ----scad-fit-----------------------------------------------------------------
fit_scad <- scad(X, y, eta = 4, minratio = 1e-2)
fit_scad

## ----mcp-fit------------------------------------------------------------------
fit_mcp <- mcp(X, y, eta = 3, minratio = 1e-2)
fit_mcp

## ----enet-fit-----------------------------------------------------------------
fit_enet <- elastic_net(X, y, lambda2 = 1, minratio = 1e-2)
fit_enet

## ----path-plots, fig.show='hold', out.width='50%'-----------------------------
fit_lasso$plot_path(xvar="fraction", log_scale = FALSE)
fit_mcp$plot_path(labels = colnames(X))
fit_scad$plot_path()
fit_enet$plot_path(standardize=FALSE)

## ----criteria-----------------------------------------------------------------
fit_lasso$criteria()

## ----criteria-plots, fig.show='hold', out.width='50%'-------------------------
fit_lasso$plot(type = "criteria")

## ----criteria-plots-2, fig.show='hold', out.width='50%'-----------------------
fit_lasso$information_criteria$plot(c("AIC", "BIC", "mBIC"))
fit_mcp$information_criteria$plot("GCV")

## ----get-model----------------------------------------------------------------
coef_lasso_bic <- fit_lasso$get_model("BIC")
coef_mcp_bic   <- fit_mcp$get_model("BIC")

cat("Non-zero Lasso coefficients (BIC):\n")
print(coef_lasso_bic[coef_lasso_bic != 0])

cat("\nNon-zero MCP coefficients (BIC):\n")
print(coef_mcp_bic[coef_mcp_bic != 0])

## ----get-model-lambda---------------------------------------------------------
lambda_lasso <- fit_lasso$major_tuning
coef_lasso_lambda <- fit_lasso$get_model(lambda_lasso[18])

cat("Non-zero Lasso coefficients for lambda = ", round(lambda_lasso[18], 2), "\n")
print(coef_lasso_lambda[coef_lasso_lambda != 0])

## ----cv, message=FALSE--------------------------------------------------------
set.seed(42)
fit_lasso$cross_validate(K = 10, verbose = FALSE)
fit_mcp$cross_validate(K = 10, verbose = FALSE)

## ----cv-plots, fig.show='hold', out.width='50%'-------------------------------
fit_lasso$plot(type = "crossval")
fit_mcp$plot(type = "crossval")

## ----cv-model-----------------------------------------------------------------
coef_lasso_cv <- fit_lasso$get_model("CV_min")
cat("Non-zero Lasso coefficients (CV min):\n")
print(coef_lasso_cv[coef_lasso_cv != 0])

## ----cv2D, message=FALSE------------------------------------------------------
set.seed(42)
fit_enet$cross_validate(K = 5, lambda2 = 10^seq(1, -3, len=20), verbose = FALSE)

## ----cv-plots-2D, fig.show='hold', out.width='50%'----------------------------
fit_enet$plot(type = "crossval")
fit_enet$cross_validation$plotCV_1D(se =  FALSE) + ggplot2::ylim(c(0.4, 0.55))

## ----predict------------------------------------------------------------------
## Predictions at the CV_min-selected model
y_hat_lasso <- fit_lasso$predict(selection = "CV_min")
y_hat_mcp   <- fit_mcp$predict(selection = "CV_min")
y_hat_enet  <- fit_enet$predict(selection = "CV_min")

## R² at the selected model
r2_lasso <- fit_lasso$r_squared[fit_lasso$get_model("CV_min", type = "index")]
r2_mcp   <- fit_mcp$r_squared[fit_mcp$get_model("CV_min", type = "index")]
r2_enet  <- fit_enet$r_squared[fit_enet$get_model("CV_min", type = "index")]
cat(sprintf("R² Lasso (CV_min): %.3f\nR² MCP   (CV_min): %.3f\nR² Enet  (CV_min): %.3f\n", r2_lasso, r2_mcp, r2_enet))

## ----stability, message=FALSE-------------------------------------------------
set.seed(42)
fit_lasso$stability(n_subsamples = 200, verbose = FALSE)

## ----stability-plot-----------------------------------------------------------
fit_lasso$plot(type = "stability")
fit_lasso$stability_path$plot(nvarsel = 7)
colnames(X)[fit_lasso$stability_path$selection(nvarsel = 7) ]

## ----debias-------------------------------------------------------------------
fit_lasso$debias <- TRUE
fit_lasso$plot_path()

## ----debias-cv----------------------------------------------------------------
fit_lasso$cross_validate(verbose = FALSE)
fit_lasso$plot(type = "crossval")
fit_lasso$debias <- FALSE
fit_lasso$plot_path()

