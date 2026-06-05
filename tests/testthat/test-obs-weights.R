context("Observation weights: consistency between quadrupen and glmnet")

tol <- 1e-2

## Shared data generated once for all tests
set.seed(2025)
n <- 100L
p  <- 50L
x <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)

## Exponential weights (sum != n, values vary widely)
w <- rexp(n, rate = 1)

test_that("unit weights give the same result as no weights", {

  fit_default <- lasso(x, y)
  fit_unit_w  <- lasso(x, y, weights = rep(1, n))

  expect_equal(as.matrix(fit_default$coefficients),
               as.matrix(fit_unit_w$coefficients),
               check.attributes = FALSE, tolerance = tol)

  expect_equal(fit_default$intercept,
               fit_unit_w$intercept,
               check.attributes = FALSE, tolerance = tol)
})

test_that("weighted lasso matches glmnet (intercept, normalized)", {

  require(glmnet)

  ## lambda1_quad = lambda_glmnet * sqrt(sum(w)):
  ## quadrupen scales by L2 norm of each column (= sqrt(n_w) * glmnet sd_w),
  ## while glmnet normalizes to unit weighted variance. The ratio gives sqrt(n_w).
  n_w      <- sum(w)
  fit_glmn <- glmnet(x, y, weights = w, lambda.min.ratio = 1e-2)
  fit_quad <- lasso(x, y, weights = w, lambda1 = fit_glmn$lambda * sqrt(n_w))

  expect_equal(as.matrix(fit_quad$coefficients),
               as.matrix(fit_glmn$beta),
               check.attributes = FALSE, tolerance = tol)

  expect_equal(fit_quad$intercept,
               fit_glmn$a0,
               check.attributes = FALSE, tolerance = tol)

  expect_equal(as.matrix(fit_quad$fitted),
               predict(fit_glmn, x),
               check.attributes = FALSE, tolerance = tol)
})

test_that("weighted lasso matches glmnet (no intercept, normalized)", {

  require(glmnet)

  n_w      <- sum(w)
  fit_glmn <- glmnet(x, y, weights = w, intercept = FALSE, lambda.min.ratio = 1e-2)
  fit_quad <- lasso(x, y, weights = w, intercept = FALSE, lambda1 = fit_glmn$lambda * sqrt(n_w))

  expect_equal(as.matrix(fit_quad$coefficients),
               as.matrix(fit_glmn$beta),
               check.attributes = FALSE, tolerance = tol)

  expect_equal(as.matrix(fit_quad$fitted),
               predict(fit_glmn, x),
               check.attributes = FALSE, tolerance = tol)
})

test_that("upweighted observations are equivalent to replicated observations", {

  ## Conceptual WLS check: upweighting observation i by integer k is exactly
  ## equivalent to having k copies of that row in the data (unit weights).
  ## Both objectives are identical term-by-term, so solutions must match
  ## at the same lambda1 values.

  upweight_idx <- 1:10
  k <- 3L

  w_up <- rep(1L, n)
  w_up[upweight_idx] <- k

  x_rep <- rbind(x, x[rep(upweight_idx, k - 1L), ])
  y_rep <- c(y, y[rep(upweight_idx, k - 1L)])

  ## Use normalize=FALSE so that column norms are the same for both fits
  ## (the replicated dataset has a different n, hence different norms otherwise).
  fit_rep <- lasso(x_rep, y_rep, normalize = FALSE, intercept = FALSE)
  fit_w   <- lasso(x, y, weights = w_up, normalize = FALSE, intercept = FALSE,
                   lambda1 = fit_rep$major_tuning)

  expect_equal(as.matrix(fit_w$coefficients),
               as.matrix(fit_rep$coefficients),
               check.attributes = FALSE, tolerance = tol)
})
