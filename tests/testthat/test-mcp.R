context("Consistency of MCP solution path (package 'ncvreg')")

testDataEnet <- readRDS("dataTest-Enet.rds")

tol <- 1e-2

get_mcp <- function(x, y, method = "quadra", plot = FALSE) {

  require(ncvreg)

  mcp_ncvreg <- ncvreg(x, y, family="gaussian", penalty = "MCP", eps = 1e-5, max.iter = 1e6)
  coef   <- mcp_ncvreg$beta[-1,]
  lambda <- mcp_ncvreg$lambda
  n <- nrow(x)

  suppressWarnings(mcp_quadru <- mcp(x, y, lambda1 = lambda * sqrt(n),
                           control = list(method = method)))
  coef_quad <- as.matrix(mcp_quadru$coefficients)

  if (plot) {
    matplot(lambda, t(coef), type = "l", log = "x")
    matplot(lambda * sqrt(n), t(coef_quad), type = "l", log = "x")
  }

  res <- list(quadru=coef_quad, ncvreg=coef)
  res
}

test_that("MCP with Newton solver", {

  ## PROSTATE DATA SET
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate

  ## Run the tests...
  out <- get_mcp(x,y)
  expect_equal(out$quadru, out$ncvreg,
               check.attributes = FALSE, tolerance = tol)

})
