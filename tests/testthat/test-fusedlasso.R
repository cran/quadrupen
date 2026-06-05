context("Consistency of the fused-Lasso")

load("dataTest-fusedlasso.rda")

tol <- 2e-4

test_that("Fused-Lasso - p=20, n=10, no intercept, Gaussian data, chain graph", {

  res <- fused_lasso(
    testData$Exn10p20X,
    testData$Exn10p20y,
    lambda1 = testData$Exn10p20Lambda1 * nrow(testData$Exn10p20X),
    lambda2 = testData$Exn10p20Lambda2 * nrow(testData$Exn10p20X), 
    intercept = FALSE,
    normalize = FALSE
  )
  expect_equal(as.vector(res$coefficients), testData$Exn10p20Sol, tolerance = tol)

})
  
test_that("Fused-Lasso - p=1000, n=20, no intercept, Gaussian data, chain graph", {

  res <- fused_lasso(
    testData$Exn20p1000X, 
    testData$Exn20p1000y,
    lambda1=testData$Exn20p1000Lambda1 * nrow(testData$Exn20p1000X),
    lambda2=testData$Exn20p1000Lambda2 * nrow(testData$Exn20p1000X), 
    intercept = FALSE,
    normalize = FALSE
  )
  sol <- testData$Exn20p1000Sol
  sol[is.na(sol)] <- res$coefficients[is.na(sol)]
  expect_equal(as.vector(res$coefficients), sol, tolerance = 1e-2)

})

test_that("Fused-Lasso - p=100, n=10, no intercept, Gaussian data, chain graph", {
  res <- fused_lasso(
    testData$Exn20p100X, 
    testData$Exn20p100y,
    lambda1 = testData$Exn20p100Lambda1 * nrow(testData$Exn20p100X), 
    lambda2 = testData$Exn20p100Lambda2 * nrow(testData$Exn20p100X), 
    intercept = FALSE,
    normalize = FALSE
  )
  sol <- testData$Exn20p100Sol
  sol[is.na(sol)] <- res$coefficients[is.na(sol)]
  expect_equal(as.vector(res$coefficients), sol, tolerance = 1e-2)
  
})

# test_that("Fused-Lasso - p=100, n=10, no intercept, Gaussian data, 2D graph", {
#   
#   G <- igraph::make_lattice(c(10,10)) |> igraph::as_adjacency_matrix(sparse = FALSE)
#   
#   res <- fused_lasso(
#     testData$Exn20p100X, 
#     testData$Exn20p100y,
#     lambda1 = testData$Exn20p100Lambda1 * nrow(testData$Exn20p100X),
#     lambda2 = testData$Exn20p100Lambda2 * nrow(testData$Exn20p100X),
#     struct = G ,
#     intercept = FALSE,
#     normalize = TRUE
#   )
#   sol <- testData$Exn20p100Sol2Dim
#   sol[is.na(sol)] <- res$coefficients[is.na(sol)]
#   expect_equal(as.vector(res$coefficients), sol, tolerance = 1e-2)
#   
# })

test_that("Fused-Lasso - p=50, n=20, no intercept, Gaussian data, chain graph", {
  res <- 
    fused_lasso(
      testData$Exn50p20X,
      testData$Exn50p20y,
      lambda1 = testData$Exn50p20Lambda1 * nrow(testData$Exn50p20X), 
      lambda2 = testData$Exn50p20Lambda2 * nrow(testData$Exn50p20X),
      intercept = FALSE,
      normalize = FALSE
    )
  sol <- testData$Exn50p20Sol
  sol[is.na(sol)] <- res$coefficients[is.na(sol)]
  expect_equal(as.vector(res$coefficients), sol, tolerance = 1e-2)

})
