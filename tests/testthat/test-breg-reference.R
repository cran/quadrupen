context("Consistency of Bounded regression")

testDataBreg <- readRDS("dataTest-boundedReg.rds")

tol <- 1e-3

test_that("Bounded regression with lambda2 = 0, QUADRA - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 0,
    control = list(method = "quadra")
    )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_0$coefficients, tolerance = tol)

})

test_that("Bounded regression with lambda2 = 5, QUADRA - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 5,
    control = list(method = "quadra")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_5$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 10 + S, QUADRA - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 10, 
    struct = testDataBreg$S,
    control = list(method = "quadra")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_10_S$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 0, QUADRA with conjuguate gradient - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 0,
    control = list(usechol = FALSE, method = "quadra")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_0$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 5, QUADRA with conjuguate gradient - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 5,
    control = list(usechol = FALSE, method = "quadra")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_5$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 10 + S, QUADRA with conjuguate gradient - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 10, 
    struct = testDataBreg$S,
    control = list(usechol = FALSE, method = "quadra")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_10_S$coefficients, tolerance = tol)
  
})


tol <- 1e-2

test_that("Bounded regression with lambda2 = 0, FISTA - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 0, control = list(method = "fista")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_0$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 5, FISTA - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 5, control = list(method = "fista")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_5$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 10 + S, FISTA  - test on the documentation example", {
  
  res <- bounded_reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 10, 
    struct = testDataBreg$S, control = list(method = "fista")
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_10_S$coefficients, tolerance = tol)
  
})

