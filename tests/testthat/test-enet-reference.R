context("Consistency of the Elastic-net solution path (package 'elasticnet')")

testDataEnet <- readRDS("dataTest-Enet.rds")

require(elasticnet)

get.enet <- function(x,y,intercept,normalize=TRUE,method="quadra",verbose=0) {
  lambda2 <- runif(1,0,10)
  enet.larsen <- enet(x,y,lambda=lambda2,intercept=intercept,normalize=normalize)
  iols <- length(enet.larsen$penalty)
  lambda1 <- enet.larsen$penalty[-iols]/2
  
  enet.quadru <- elastic_net(x,y,intercept=intercept,normalize=normalize,
                             lambda1=lambda1, lambda2=lambda2,
                             control = list(method=method,verbose=verbose))
  
  quad <- as.matrix(enet.quadru$coefficients)
  
  enet <- predict(enet.larsen, type="coefficients", naive=TRUE)$coefficients[-iols,]
  
  return(list(quad = quad, enet = t(enet)))
}

test_that("Elastic-net is correct w.r.t a reference solution", {

  tol <- 1e-3

  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate
  
  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## RANDOM DATA
  x <- testDataEnet$x_sim
  y <- testDataEnet$y_sim
  
  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

})

test_that("Elastic-net is correct w.r.t a reference solution - FISTA", {

  tol <- 1e-2
  
  ## PROSTATE DATA SET
  x <- testDataEnet$x_sim[, 1:20]
  y <- testDataEnet$y_sim
  
  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE, method="fista")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)
  
  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE, method="fista")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## Run the tests...
  without.intercept <- get.enet(x,y,intercept=FALSE, method="fista")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)
  
  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE, method="fista")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)
  

})

test_that("Elastic-net is correct w.r.t a reference solution - PGD", {

  tol <- 1e-2

  ## PROSTATE DATA SET
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate

  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE, method="pgd")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE, method="pgd")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## Run the tests...
  without.intercept <- get.enet(x,y,intercept=FALSE, method="pgd")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE, method="pgd")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)

})
