context("Consistency of the Group-Lasso solution path")

testData <- readRDS("dataTest-GroupLasso.rds")

tol <- 1e-2

get_grplasso <- function(x, y, group, lambda, intercept, normalize, method = "quadra")  {
  
  out_quadr   <- quadrupen::group_lasso(x, y, group, lambda1 = lambda, 
                                        intercept = intercept, normalize = normalize,
                                        control = list(method = method))
  coef_quadr  <-  as.matrix(out_quadr$coefficients)
  group_quadr <- rowsum(coef_quadr^2, group)
  inter_quadr <- out_quadr$intercept
  
  res <- list(coef = coef_quadr, group = group_quadr, intercept = inter_quadr, lambda = lambda)
  res
}

test_that("Group-Lasso with lambda2 = 0, intercept and normalization, QUADRA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_inter_norm$lambda, TRUE, TRUE)
  
  expect_equal(quad, testData$grplasso_inter_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, intercept and no normalization - QUADRA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_inter_nonorm$lambda, TRUE, FALSE)
  
  expect_equal(quad, testData$grplasso_inter_nonorm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, no intercept and normalization - QUADRA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_norm$lambda, FALSE, TRUE)
  
  expect_equal(quad, testData$grplasso_nointer_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, no intercept and no normalization - QUADRA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_nonorm$lambda, FALSE, FALSE)
  
  expect_equal(quad, testData$grplasso_nointer_nonorm, check.attributes = FALSE, tolerance = tol)
}
)



test_that("Group-Lasso with lambda2 = 0, intercept and normalization, FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_inter_norm$lambda, TRUE, TRUE, "fista")
  
  expect_equal(quad, testData$grplasso_inter_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, intercept and no normalization - FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_inter_nonorm$lambda, TRUE, FALSE, "fista")
  
  expect_equal(quad, testData$grplasso_inter_nonorm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, no intercept and normalization - FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_norm$lambda, FALSE, TRUE, "fista")
  
  expect_equal(quad, testData$grplasso_nointer_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, no intercept and no normalization - FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_nonorm$lambda, FALSE, FALSE, "fista")
  
  expect_equal(quad, testData$grplasso_nointer_nonorm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Quadratic solver is facter than FISTA", {

  time_quadra <- system.time(
    quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_nonorm$lambda, FALSE, FALSE, "quadra")  
  )[[3]]
  
  time_fista <- system.time(
    quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_nonorm$lambda, FALSE, FALSE, "fista")  
  )[[3]]
  
  expect_lt(time_quadra, time_fista)
  }
)

