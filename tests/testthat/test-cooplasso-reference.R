context("Consistency of the Cooperative-Lasso solution path")

testData <- readRDS("dataTest-CoopLasso.rds")

tol <- 1e-2

get_cooplasso <- function(x, y, group, lambda, intercept, normalize, method = "fista") {
  
  out_quadr   <- quadrupen::coop_lasso(x, y, group, lambda1 = lambda, 
                                      intercept = intercept, normalize = normalize,
                                      control = list(method = method))
  coef_quadr  <-  as.matrix(out_quadr$coefficients)
  group_quadr <- rowsum(coef_quadr^2, group)
  inter_quadr <- out_quadr$intercept
  
  res <- list(coef = coef_quadr, group = group_quadr, intercept = inter_quadr, lambda = lambda)
  res
}

test_that("Coop-Lasso with lambda2 = 0, intercept and normalization, FISTA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_inter_norm$lambda, TRUE, TRUE)
  
  expect_equal(quad, testData$cooplasso_inter_norm, check.attributes = FALSE, tolerance = tol)
  }
)

test_that("Coop-Lasso with lambda2 = 0, intercept and no normalization - FISTA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_inter_nonorm$lambda, TRUE, FALSE)
  
  expect_equal(quad, testData$cooplasso_inter_nonorm, check.attributes = FALSE, tolerance = tol)
  }
)

test_that("Coop-Lasso with lambda2 = 0, no intercept and normalization - FISTA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_nointer_norm$lambda, FALSE, TRUE)
  
  expect_equal(quad, testData$cooplasso_nointer_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Coop-Lasso with lambda2 = 0, no intercept and no normalization - FISTA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_nointer_nonorm$lambda, FALSE, FALSE)
  
  expect_equal(quad, testData$cooplasso_nointer_nonorm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Coop-Lasso with lambda2 = 0, intercept and normalization, QUADRA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_inter_norm$lambda, TRUE, TRUE, "quadra")
  
  expect_equal(quad, testData$cooplasso_inter_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Coop-Lasso with lambda2 = 0, intercept and no normalization - QUADRA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_inter_nonorm$lambda, TRUE, FALSE, "quadra")
  
  expect_equal(quad, testData$cooplasso_inter_nonorm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Coop-Lasso with lambda2 = 0, no intercept and normalization - QUADRA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_nointer_norm$lambda, FALSE, TRUE, "quadra")
  
  expect_equal(quad, testData$cooplasso_nointer_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Coop-Lasso with lambda2 = 0, no intercept and no normalization - QUADRA - test on the documentation example", {
  quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$cooplasso_nointer_nonorm$lambda, FALSE, FALSE, "quadra")
  
  expect_equal(quad, testData$cooplasso_nointer_nonorm, check.attributes = FALSE, tolerance = tol)
}
)

# test_that("Quadratic solver is facter than FISTA", {
# 
#   time_quadra <- system.time(
#     quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_nonorm$lambda, FALSE, FALSE, "quadra")
#   )[[3]]
# 
#   time_fista <- system.time(
#     quad <- get_cooplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_nonorm$lambda, FALSE, FALSE, "fista")
#   )[[3]]
# 
#   expect_lt(time_quadra, time_fista)
# }
# )

