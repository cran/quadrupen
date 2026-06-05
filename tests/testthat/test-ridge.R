context("Consistency of the Ridge solution path (package MASS)")

testDataEnet <- readRDS("dataTest-Enet.rds")

tol <- 1e-4

require(MASS)

test_that("lasso_quad2lars", {


  get_ridge <- function(x,y,intercept) {

    lambda <- 10^seq(2,-2,len = 100)
    if (intercept) {
      ridge_mass <- MASS::lm.ridge(y ~ x + 1, lambda = lambda * nrow(x))
      mass <- list(coef = t(coefficients(ridge_mass)[, -1]),
                   mu   = coefficients(ridge_mass)[, 1])
    } else {
      ridge_mass <- MASS::lm.ridge(y ~ x + 0, lambda = lambda * nrow(x))
      mass <- list(coef = t(coefficients(ridge_mass)),
                   mu   = rep(0,length(lambda)))

    }
    ridge_quad <- ridge(x, y, intercept=intercept, normalize=TRUE, lambda=lambda)
    quad <- list(coef = as.matrix(ridge_quad$coefficients),
                 mu   = ridge_quad$intercept)


    return(list(quad=quad,mass=mass))
  }

  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate

  ## Run the tests...
  with.intercept <- get_ridge(x,y,TRUE)
  expect_equal(with.intercept$quad,
               with.intercept$mass, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get_ridge(x,y,FALSE)
  expect_equal(without.intercept$quad,
               without.intercept$mass, check.attributes = FALSE, tolerance = tol)

  ## RANDOM DATA
  x <- testDataEnet$x_sim
  y <- testDataEnet$y_sim

  ## Run the tests...
  ## Run the tests...
  with.intercept <- get_ridge(x,y,TRUE)
  expect_equal(with.intercept$quad,
               with.intercept$mass, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get_ridge(x,y,FALSE)
  expect_equal(without.intercept$quad,
               without.intercept$mass, check.attributes = FALSE, tolerance = tol)
})
