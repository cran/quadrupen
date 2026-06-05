context("Independent verification of SCAD and MCP solvers")

# ─────────────────────────────────────────────────────────────────────────────
# Pure-R reference implementations
# ─────────────────────────────────────────────────────────────────────────────

# SCAD proximal operator (scalar then vector via Vectorize)
.scad_prox <- function(z, lambda, a) {
  abs_z <- abs(z)
  dplyr_free <- function(zi) {
    ai <- abs(zi)
    if      (ai <= lambda)          0
    else if (ai <= 2 * lambda)      sign(zi) * (ai - lambda)
    else if (ai <= a   * lambda)    sign(zi) * ((a - 1) * ai - a * lambda) / (a - 2)
    else                            zi
  }
  vapply(z, dplyr_free, numeric(1))
}

# MCP proximal operator
.mcp_prox <- function(z, lambda, gamma) {
  vapply(z, function(zi) {
    ai <- abs(zi)
    if      (ai <= lambda)          0
    else if (ai <= gamma * lambda)  sign(zi) * gamma * (ai - lambda) / (gamma - 1)
    else                            zi
  }, numeric(1))
}

# SCAD penalty value (sum over coordinates)
.scad_pen <- function(beta, lambda, a) {
  ab <- abs(beta)
  z1 <- ab <= lambda
  z2 <- !z1 & ab <= a * lambda
  z3 <- ab  > a * lambda
  pen <- numeric(length(beta))
  pen[z1] <- lambda * ab[z1]
  pen[z2] <- -(ab[z2]^2 - 2 * a * lambda * ab[z2] + lambda^2) / (2 * (a - 1))
  pen[z3] <- (a + 1) * lambda^2 / 2
  sum(pen)
}

# MCP penalty value
.mcp_pen <- function(beta, lambda, gamma) {
  ab <- abs(beta)
  z1 <- ab <= gamma * lambda
  z2 <- ab  > gamma * lambda
  pen <- numeric(length(beta))
  pen[z1] <- lambda * ab[z1] - ab[z1]^2 / (2 * gamma)
  pen[z2] <- gamma * lambda^2 / 2
  sum(pen)
}

# Criterion 1/2 ||y - Xβ||² + λ·pen(β)  (quadrupen scale, no intercept)
.criterion <- function(beta, x, y, lambda, eta, type) {
  loss <- 0.5 * sum((y - x %*% beta)^2)
  pen  <- if (type == "scad") .scad_pen(beta, lambda, eta) else .mcp_pen(beta, lambda, eta)
  loss + pen
}

# KKT check for normalize=FALSE, intercept=FALSE, lambda2=0
# Returns max KKT violation (should be ≈ 0 at a true stationary point).
.kkt_violation <- function(beta, x, y, lambda, eta, type) {
  grad <- drop(t(x) %*% (x %*% beta - y))
  ab   <- abs(beta)
  inactive <- ab < 1e-8

  # Inactive KKT: |grad_j| ≤ λ
  viol_inact <- if (any(inactive))  max(abs(grad[inactive]) - lambda, 0) else 0

  # Active KKT: |grad_j + p'(|β_j|)·sign(β_j)| = 0
  if (!any(!inactive)) return(viol_inact)

  bA <- beta[!inactive]; gA <- grad[!inactive]; abA <- ab[!inactive]
  if (type == "scad") {
    d <- ifelse(abA <= lambda,            lambda,
         ifelse(abA <= eta * lambda,      (eta * lambda - abA) / (eta - 1),
                                          0))
  } else { # MCP
    d <- ifelse(abA <= eta * lambda,      lambda - abA / eta,
                                          0)
  }
  viol_act <- max(abs(gA + d * sign(bA)))
  max(viol_inact, viol_act)
}

# ─────────────────────────────────────────────────────────────────────────────
# Test 1 — Orthogonal design: solution == proximal operator applied to X^T y
# ─────────────────────────────────────────────────────────────────────────────
# When X has orthonormal columns (X^T X = I_p) and normalize = FALSE,
# intercept = FALSE, lambda2 = 0, the global minimiser is prox_λ(X^T y)
# applied coordinate-by-coordinate. No external package needed.

tol_orth <- 1e-4

test_that("SCAD — orthogonal design matches analytical proximal", {
  set.seed(2024)
  n <- 80 ; p <- 30
  X <- qr.Q(qr(matrix(rnorm(n * p), n, p)))  # X^T X = I_p exactly
  beta_true <- c(0, 3, -2, 0, 5, rep(0, p - 5))
  y <- X %*% beta_true + rnorm(n, 0, 0.3)

  lambda <- 0.4 ; eta <- 3.7
  z <- drop(t(X) %*% y)   # OLS == X^T y when X^T X = I_p
  beta_ref <- .scad_prox(z, lambda, eta)

  suppressWarnings(
    fit <- scad(X, y, lambda1 = lambda, lambda2 = 0,
                intercept = FALSE, normalize = FALSE,
                control = list(threshold = 1e-9, maxiter = 500))
  )
  beta_quad <- drop(as.matrix(fit$coefficients))

  expect_equal(beta_quad, beta_ref, tolerance = tol_orth, check.attributes = FALSE)
})

test_that("MCP — orthogonal design matches analytical proximal", {
  set.seed(2024)
  n <- 80 ; p <- 30
  X <- qr.Q(qr(matrix(rnorm(n * p), n, p)))
  beta_true <- c(0, 3, -2, 0, 5, rep(0, p - 5))
  y <- X %*% beta_true + rnorm(n, 0, 0.3)

  lambda <- 0.4 ; eta <- 3
  z <- drop(t(X) %*% y)
  beta_ref <- .mcp_prox(z, lambda, eta)

  suppressWarnings(
    fit <- mcp(X, y, lambda1 = lambda, lambda2 = 0,
               intercept = FALSE, normalize = FALSE,
               control = list(threshold = 1e-9, maxiter = 500))
  )
  beta_quad <- drop(as.matrix(fit$coefficients))

  expect_equal(beta_quad, beta_ref, tolerance = tol_orth, check.attributes = FALSE)
})

# ─────────────────────────────────────────────────────────────────────────────
# Test 2 — KKT conditions: every solution on the path must be a stationary point
# ─────────────────────────────────────────────────────────────────────────────
# This test is package-independent: it only uses X, y, and the definition of
# the SCAD/MCP subdifferential. Works for any dimension.

tol_kkt <- 1e-3   # should be near machine precision × cond(X); relax for poorly scaled problems

run_kkt_test <- function(x, y, type, eta, tol) {
  fn <- if (type == "scad") scad else mcp
  suppressWarnings(
    fit <- fn(x, y, eta = eta, lambda2 = 0,
              intercept = FALSE, normalize = FALSE,
              control = list(threshold = 1e-8, maxiter = 500))
  )
  lambdas <- fit$major_tuning
  coefs   <- as.matrix(fit$coefficients)

  viols <- vapply(seq_along(lambdas), function(k)
    .kkt_violation(coefs[, k], x, y, lambdas[k], eta, type),
    numeric(1))

  max(viols)
}

test_that("SCAD — KKT conditions satisfied on entire path (small data)", {
  testData <- readRDS("dataTest-Enet.rds")
  x <- testData$x_prostate ; y <- testData$y_prostate
  viol <- run_kkt_test(x, y, "scad", eta = 3.7, tol = tol_kkt)
  expect_lt(viol, tol_kkt)
})

test_that("MCP — KKT conditions satisfied on entire path (small data)", {
  testData <- readRDS("dataTest-Enet.rds")
  x <- testData$x_prostate ; y <- testData$y_prostate
  viol <- run_kkt_test(x, y, "mcp", eta = 3, tol = tol_kkt)
  expect_lt(viol, tol_kkt)
})

test_that("SCAD — KKT conditions satisfied on entire path (larger data)", {
  testData <- readRDS("dataTest-Enet.rds")
  x <- testData$x_sim ; y <- testData$y_sim
  viol <- run_kkt_test(x, y, "scad", eta = 3.7, tol = tol_kkt)
  expect_lt(viol, tol_kkt)
})

test_that("MCP — KKT conditions satisfied on entire path (larger data)", {
  testData <- readRDS("dataTest-Enet.rds")
  x <- testData$x_sim ; y <- testData$y_sim
  viol <- run_kkt_test(x, y, "mcp", eta = 3, tol = tol_kkt)
  expect_lt(viol, tol_kkt)
})

# ─────────────────────────────────────────────────────────────────────────────
# Test 3 — Objective value: quadrupen ≤ ncvreg (when they disagree, ours wins)
# ─────────────────────────────────────────────────────────────────────────────
# Both methods should find a (local) minimum. If they return different
# solutions, the one with the lower criterion value is preferable.
# Uses normalize = FALSE so objective scales are directly comparable.
#
# ncvreg criterion: (1/2n)||y-Xβ||² + λ pen(β)  → multiply by n to match quadrupen
# quadrupen uses lambda1 = lambda_ncvreg * n  (without normalisation factor)

test_that("SCAD — objective value not worse than ncvreg (larger data)", {
  require(ncvreg)
  testData <- readRDS("dataTest-Enet.rds")
  x <- testData$x_sim ; y <- testData$y_sim
  n <- nrow(x)

  # Fit ncvreg (without standardisation to keep scales matching)
  suppressWarnings(fit_nc <- ncvreg(x, y, family = "gaussian", penalty = "SCAD",
                   eps = 1e-8, max.iter = 1e6))
  lambdas_nc <- fit_nc$lambda

  # Fit quadrupen at the same effective lambda (multiply by n to match quadrupen scale)
  suppressWarnings(
    fit_qu <- scad(x, y, lambda1 = lambdas_nc * n, lambda2 = 0,
                   intercept = FALSE, normalize = FALSE,
                   control = list(threshold = 1e-8, maxiter = 500))
  )

  coef_nc <- fit_nc$beta[-1, ]  # drop intercept row
  coef_qu <- as.matrix(fit_qu$coefficients)

  eta <- 3.7
  # Compare at each lambda: quadrupen objective ≤ ncvreg objective + small slack
  # (both evaluated on the quadrupen-scale criterion)
  diffs <- vapply(seq_along(lambdas_nc), function(k) {
    lam <- lambdas_nc[k] * n
    obj_nc <- .criterion(coef_nc[, k], x, y, lam, eta, "scad")
    obj_qu <- .criterion(coef_qu[, k], x, y, lam, eta, "scad")
    obj_qu - obj_nc   # negative means quadrupen is better
  }, numeric(1))

  # Quadrupen should never be worse (up to numerical tolerance)
  expect_true(all(diffs <= 1e-4),
    info = paste0("quadrupen worse than ncvreg at ", sum(diffs > 1e-4),
                  " lambda values; max excess = ", round(max(diffs), 6)))
})

test_that("MCP — objective value not worse than ncvreg (larger data)", {
  require(ncvreg)
  testData <- readRDS("dataTest-Enet.rds")
  x <- testData$x_sim ; y <- testData$y_sim
  n <- nrow(x)

  suppressWarnings(
  fit_nc <- ncvreg(x, y, family = "gaussian", penalty = "MCP",
                   eps = 1e-8, max.iter = 1e6)
  )
  lambdas_nc <- fit_nc$lambda

  suppressWarnings(
    fit_qu <- mcp(x, y, lambda1 = lambdas_nc * n, lambda2 = 0,
                  intercept = FALSE, normalize = FALSE,
                  control = list(threshold = 1e-8, maxiter = 500))
  )

  coef_nc <- fit_nc$beta[-1, ]
  coef_qu <- as.matrix(fit_qu$coefficients)

  eta <- 3
  diffs <- vapply(seq_along(lambdas_nc), function(k) {
    lam <- lambdas_nc[k] * n
    obj_nc <- .criterion(coef_nc[, k], x, y, lam, eta, "mcp")
    obj_qu <- .criterion(coef_qu[, k], x, y, lam, eta, "mcp")
    obj_qu - obj_nc
  }, numeric(1))

  expect_true(all(diffs <= 1e-4),
    info = paste0("quadrupen worse than ncvreg at ", sum(diffs > 1e-4),
                  " lambda values; max excess = ", round(max(diffs), 6)))
})
