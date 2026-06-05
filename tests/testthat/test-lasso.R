context("Consistency of the Lasso solution paths (package 'lars' and 'glmnet')")

testDataEnet <- readRDS("dataTest-Enet.rds")

tol <- 1e-2

test_that("lasso_quad2lars", {

  require(lars)

  get.lars <- function(x,y,intercept,normalize) {
      lasso.larsen <- lars(x,y,intercept=intercept,normalize=normalize)
      iols <- nrow(lasso.larsen$beta) ## remove last entry corresponding to the OLS estimator
      lambda1 <-  lasso.larsen$lambda ## usde the lars lambda grid
      lasso.quadru <- elastic_net(x,y, intercept=intercept, normalize=normalize,
                                  lambda1=lambda1, lambda2=0, control=list(method="quadra"))
      quad <- list(coef = as.matrix(lasso.quadru$coefficients),
                   rss  = deviance(lasso.quadru))

      lars <- list(coef = t(lasso.larsen$beta[-iols, ]),
                   rss  = lasso.larsen$RSS[-iols])

      return(list(quad=quad,lars=lars))
  }

  ## PROSTATE DATA SET
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate

  ## Run the tests...
  with.intercept <- get.lars(x,y,TRUE,TRUE)
  expect_equal(with.intercept$quad,
               with.intercept$lars, check.attributes = FALSE, tolerance = tol)

  with.intercept.unnormalized <-get.lars(x,y,TRUE,FALSE)
  expect_equal(with.intercept.unnormalized$quad,
               with.intercept.unnormalized$lars, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.lars(x,y,FALSE,TRUE)
  expect_equal(without.intercept$quad,
               without.intercept$lars, check.attributes = FALSE, tolerance = tol)

  without.intercept.unnormalized <- get.lars(x,y,FALSE,FALSE)
  expect_equal(without.intercept.unnormalized$quad,
               without.intercept.unnormalized$lars, check.attributes = FALSE, tolerance = tol)

  ## RANDOM DATA
  x <- testDataEnet$x_sim
  y <- testDataEnet$y_sim

  ## Run the tests...
  with.intercept <-get.lars(x,y,TRUE,TRUE)
  expect_equal(with.intercept$coef.quad,
               with.intercept$coef.lars, check.attributes = FALSE, tolerance = tol)

  with.intercept.unnormalized <-get.lars(x,y,TRUE,FALSE)
  expect_equal(with.intercept.unnormalized$coef.quad,
               with.intercept.unnormalized$coef.lars, check.attributes = FALSE, tolerance = tol)

  without.intercept <-get.lars(x,y,FALSE,TRUE)
  expect_equal(without.intercept$coef.quad,
               without.intercept$coef.lars, check.attributes = FALSE, tolerance = tol)

  without.intercept.unnormalized <-get.lars(x,y,FALSE,FALSE)
  expect_equal(without.intercept.unnormalized$coef.quad,
               without.intercept.unnormalized$coef.lars, check.attributes = FALSE, tolerance = tol)

})

test_that("lasso_quad2glmnet", {

  require(glmnet)

  ## SECOND CHECK: compare to glmnet with prescaling of x
  x <- matrix(rnorm(100*50),100,50)
  y <- rnorm(100)
  y <- y - mean(y)
  n <- nrow(x)
  p <- ncol(x)

  ## If thresh is set to the default, the test won't pass!!!
  ## This is because coordinate descent is fast yet not extremely accurate
  lasso.glmn <- glmnet(x,y, lambda.min.ratio=1e-2)#, thresh=1e-20)
  lasso.quad <- elastic_net(x,y, lambda1=lasso.glmn$lambda*sqrt(n), lambda2=0)

  quad <- list(coef   = as.matrix(lasso.quad$coefficients),
               fitted = as.matrix(lasso.quad$fitted))

  glmn <- list(coef   = as.matrix(lasso.glmn$beta),
               fitted = predict(lasso.glmn,x))


  expect_equal(quad, glmn, check.attributes = FALSE, tolerance = tol)

})
