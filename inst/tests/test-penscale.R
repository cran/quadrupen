context("Brief consistency of penscale")

test_that("Consistency between weighted 'quadrupen' and weighted sof-thresholding (X=Identity)", {

  ## SIMPLE CHECK: used identity for design matrix
  n <- 100
  p <- 100
  x <- diag(rep(1,n))
  y <- rnorm(100)
  ## no penscale...
  lasso.quad <- elastic.net(x,y, intercept=FALSE, lambda2=0)
  theo.path  <- t(sapply(lasso.quad@lambda1, function(lambda) y*pmax(0,1-lambda/abs(y))))

  expect_that(as.matrix(lasso.quad@coefficients), is_equivalent_to(theo.path))

  ## random penscale...
  w <- 1/runif(p,0.5,1)
  lasso.quad <- elastic.net(x,y, intercept=FALSE, penscale=w, lambda2=0)
  theo.path  <- t(sapply(lasso.quad@lambda1, function(lambda) y*pmax(0,1-lambda*w/abs(y))))

  expect_that(as.matrix(lasso.quad@coefficients), is_equivalent_to(theo.path))
})

