context("Consistency of LAVA")

testDataEnet <- readRDS("dataTest-Enet.rds")

tol <- 1e-4

test_that("lava as lasso (when lambda2 -> infty) ", {
  
  
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate
  
  ## INTERCEPT, NORMALIZE
  lasso_out <- lasso(x, y)
  lava_out  <- lava(x, y, lambda1 = lasso_out$major_tuning, lambda2 = 1E5)

  with_intercept_lasso <- list(coef = as.matrix(lasso_out$coefficients), mu = lasso_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_lasso,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)

  ## NO INTERCEPT, NORMALIZE
  lasso_out <- lasso(x, y, intercept = FALSE)
  lava_out  <- lava(x, y, intercept = FALSE, lambda1 = lasso_out$major_tuning, lambda2 = 1E5)
  
  with_intercept_lasso <- list(coef = as.matrix(lasso_out$coefficients), mu = lasso_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_lasso,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)

  ## INTERCEPT, NO NORMALIZE
  lasso_out <- lasso(x, y, normalize = FALSE)
  lava_out  <- lava(x, y, normalize = FALSE, lambda1 = lasso_out$major_tuning, lambda2 = 1E8)
  
  with_intercept_lasso <- list(coef = as.matrix(lasso_out$coefficients), mu = lasso_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_lasso,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)
  
  ## NO INTERCEPT, NO NORMALIZE
  lasso_out <- lasso(x, y, normalize = FALSE, intercept = FALSE)
  lava_out  <- lava(x, y, normalize = FALSE, intercept = FALSE, lambda1 = lasso_out$major_tuning, lambda2 = 1E8)
  
  with_intercept_lasso <- list(coef = as.matrix(lasso_out$coefficients), mu = lasso_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_lasso,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)
  
})


test_that("lava as ridge (when lambda1 -> infty) ", {
  
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate
  
  ## INTERCEPT, NORMALIZE
  ridge_out <- ridge(x, y, lambda = .75)
  lava_out  <- lava(x, y, lambda1 = 1e8, lambda2 = .75)
  with_intercept_ridge <- list(coef = as.matrix(ridge_out$coefficients), mu = ridge_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_ridge,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)
  
  ## NO INTERCEPT, NORMALIZE
  ridge_out <- ridge(x, y, intercept = FALSE, lambda = 2.75)
  lava_out  <- lava(x, y, intercept = FALSE, lambda1 = 1e8, lambda2 = 2.75)
  with_intercept_ridge <- list(coef = as.matrix(ridge_out$coefficients), mu = ridge_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_ridge,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)
  
  ## INTERCEPT, NO NORMALIZE
  ridge_out <- ridge(x, y, normalize = FALSE, lambda = 1.75)
  lava_out  <- lava(x, y, normalize = FALSE, lambda1 = 1e8, lambda2 = 1.75)
  with_intercept_ridge <- list(coef = as.matrix(ridge_out$coefficients), mu = ridge_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_ridge,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)
  
  ## NO INTERCEPT, NO NORMALIZE
  ridge_out <- ridge(x, y, intercept = FALSE, normalize = FALSE, lambda = .5)
  lava_out  <- lava(x, y, intercept = FALSE, normalize = FALSE, lambda1 = 1e8, lambda2 = .5)
  with_intercept_ridge <- list(coef = as.matrix(ridge_out$coefficients), mu = ridge_out$intercept)
  with_intercept_lava  <- list(coef = lava_out$coefficients, mu = lava_out$intercept)
  expect_equal(with_intercept_ridge,
               with_intercept_lava , check.attributes = FALSE, tolerance = tol)
  
})
