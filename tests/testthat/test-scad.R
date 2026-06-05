context("Consistency of SCAD solution path (package 'ncvreg')")

testDataEnet <- readRDS("dataTest-Enet.rds")

tol <- 1e-2

get_scad <- function(x, y, method = "quadra", plot = FALSE) {

  require(ncvreg)

  scad_ncvreg <- ncvreg(x, y, family = "gaussian", penalty = "SCAD", eps = 1e-5, max.iter = 1e6)
  coef   <- scad_ncvreg$beta[-1,]

  lambda <- scad_ncvreg$lambda ; n <- nrow(x)
  suppressWarnings(scad_quadru <- scad(x, y,  lambda1 = lambda * sqrt(n),
                      control = list(method = method)))
  coef_quad <- as.matrix(scad_quadru$coefficients)

  if (plot) {
    matplot(lambda, t(coef), type = "l", log = "x")
    matplot(lambda * sqrt(n), t(coef_quad), type = "l", log = "x")
  }

  res <- list(quadru=coef_quad, ncvreg=coef)
  res
}

test_that("SCAD with Quadratic solver", {

  ## PROSTATE DATA SET
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate

  ## Run the tests...
  out <- get_scad(x,y)
  expect_equal(out$quadru, out$ncvreg,
               check.attributes = FALSE, tolerance = tol)

})
