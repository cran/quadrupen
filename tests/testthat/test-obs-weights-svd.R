context("Observation weights for SVD-based methods (Ridge, Lava, GroupLava)")

tol_ridge <- 1e-3
tol_lava  <- 1e-2

set.seed(2025)
n <- 100L
p  <- 50L
x <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)

w <- rexp(n, rate = 1)   # general weights, sum != n

## ── Ridge ─────────────────────────────────────────────────────────────────────

test_that("weighted ridge: unit weights give the same result as no weights", {

  fit_default <- ridge(x, y)
  fit_unit_w  <- ridge(x, y, weights = rep(1, n))

  expect_equal(as.matrix(fit_default$coefficients),
               as.matrix(fit_unit_w$coefficients),
               check.attributes = FALSE, tolerance = tol_ridge)
  expect_equal(fit_default$intercept,
               fit_unit_w$intercept,
               check.attributes = FALSE, tolerance = tol_ridge)
})

test_that("weighted ridge matches direct WLS formula", {

  ## Direct WLS ridge closed form (no external package):
  ##   beta_hat(lambda) = (X_n^T W X_n + lambda I)^{-1} X_n^T W y_c / norm_X
  ## where X_n is column-standardized with weighted norms, y_c is weighted-mean centred.
  ## MASS::lm.ridge does NOT implement standard WLS ridge, so we avoid it.

  lam <- 10^seq(1, -1, length.out = 20)
  n_w    <- sum(w)
  x_bar  <- as.numeric(t(x) %*% w / n_w)
  y_bar  <- sum(w * y) / n_w
  xc     <- sweep(x, 2, x_bar)
  normx  <- sqrt(as.numeric(t(w) %*% (xc^2)))
  xn     <- sweep(xc, 2, normx, "/")
  yc     <- y - y_bar
  XtWX   <- crossprod(xn, w * xn)
  XtWy   <- as.numeric(crossprod(xn, w * yc))

  direct <- sapply(lam, function(l) {
    b_norm <- solve(XtWX + l * diag(p), XtWy)
    b_norm / normx
  })
  mu_direct <- y_bar - as.numeric(t(x_bar) %*% direct)

  fit_quad <- ridge(x, y, weights = w, lambda = lam, intercept = TRUE, normalize = TRUE)

  expect_equal(as.matrix(fit_quad$coefficients), direct,
               check.attributes = FALSE, tolerance = tol_ridge)
  expect_equal(fit_quad$intercept, mu_direct,
               check.attributes = FALSE, tolerance = tol_ridge)
})

test_that("weighted ridge: upweighted observations equivalent to replicated observations", {

  upweight_idx <- 1:10
  k <- 3L

  w_up   <- rep(1L, n) ; w_up[upweight_idx] <- k
  x_rep  <- rbind(x, x[rep(upweight_idx, k - 1L), ])
  y_rep  <- c(y, y[rep(upweight_idx, k - 1L)])

  lam <- 10^seq(1, -1, length.out = 20)
  fit_w   <- ridge(x, y,      weights = w_up, lambda = lam, normalize = FALSE)
  fit_rep <- ridge(x_rep, y_rep,              lambda = lam, normalize = FALSE)

  expect_equal(as.matrix(fit_w$coefficients),
               as.matrix(fit_rep$coefficients),
               check.attributes = FALSE, tolerance = tol_ridge)
})

## ── Lava ──────────────────────────────────────────────────────────────────────

test_that("weighted lava: unit weights give the same result as no weights", {

  fit_default <- lava(x, y)
  fit_unit_w  <- lava(x, y, weights = rep(1, n))

  expect_equal(fit_default$coefficients,
               fit_unit_w$coefficients,
               check.attributes = FALSE, tolerance = tol_lava)
})

test_that("weighted lava as weighted lasso (lambda2 -> infty)", {

  ## When the dense penalty dominates, the dense component vanishes
  ## and Lava reduces to Lasso — regardless of observation weights.
  lasso_out <- lasso(x, y, weights = w)
  lava_out  <- lava(x, y, weights = w,
                    lambda1 = lasso_out$major_tuning, lambda2 = 1e5)

  expect_equal(lava_out$coefficients,
               as.matrix(lasso_out$coefficients),
               check.attributes = FALSE, tolerance = tol_lava)
  expect_equal(lava_out$intercept,
               lasso_out$intercept,
               check.attributes = FALSE, tolerance = tol_lava)
})

test_that("weighted lava as weighted ridge (lambda1 -> infty)", {

  ## When the sparse penalty dominates, the sparse component vanishes
  ## and Lava reduces to Ridge — the most important consistency check
  ## for the weighted SVD preprocessing.
  lambda2_val <- 0.75
  ridge_out <- ridge(x, y, weights = w, lambda = lambda2_val)
  lava_out  <- lava(x, y, weights = w, lambda1 = 1e8, lambda2 = lambda2_val)

  expect_equal(lava_out$coefficients,
               as.matrix(ridge_out$coefficients),
               check.attributes = FALSE, tolerance = tol_lava)
  expect_equal(lava_out$intercept,
               ridge_out$intercept,
               check.attributes = FALSE, tolerance = tol_lava)
})

## ── GroupLava ─────────────────────────────────────────────────────────────────

test_that("weighted group_lava: unit weights give the same result as no weights", {

  grp <- rep(1:(p / 5), each = 5)

  fit_default <- group_lava(x, y, group = grp)
  fit_unit_w  <- group_lava(x, y, group = grp, weights = rep(1, n))

  expect_equal(fit_default$coefficients,
               fit_unit_w$coefficients,
               check.attributes = FALSE, tolerance = tol_lava)
})

test_that("weighted group_lava as weighted ridge (lambda1 -> infty)", {

  grp         <- rep(1:(p / 5), each = 5)
  lambda2_val <- 1.0
  ridge_out   <- ridge(x, y, weights = w, lambda = lambda2_val)
  glava_out   <- group_lava(x, y, group = grp, weights = w,
                             lambda1 = 1e8, lambda2 = lambda2_val)

  expect_equal(glava_out$coefficients,
               as.matrix(ridge_out$coefficients),
               check.attributes = FALSE, tolerance = tol_lava)
  expect_equal(glava_out$intercept,
               ridge_out$intercept,
               check.attributes = FALSE, tolerance = tol_lava)
})
