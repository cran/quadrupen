context("Consistency of the Structured Elastic-net (reference is computed via the 'augmented data' approach)")

testDataEnet <- readRDS("dataTest-Enet.rds")

tol <- 1e-3

test_that("Consistency of the structured elastic-net", {

  require(elasticnet)

  get.coef <- function(x,y,intercept=TRUE,normalize=TRUE,C=diag(rep(1,ncol(x)))) {
    lambda2 <- runif(1,0,10)

    ## INTERCEPT AND NORMALIZATION TREATMENT
    if (intercept) {
      xbar <- colMeans(x)
      ybar <- mean(y)
    } else {
      xbar <- rep(0,p)
      ybar <- 0
    }
    xs <- scale(x,xbar,FALSE)
    ys <- y-ybar

    if (normalize) {
      normx <- sqrt(drop(colSums(xs^2)))
      xs <- scale(xs,FALSE,normx)
    } else {
      normx <- rep(1,p)
    }

    ## The reference: augmented data approach
    x.aug <- rbind(xs,sqrt(lambda2)*C)
    y.aug <- c(ys,rep(0,ncol(x)))
    senet1 <- enet(x.aug,y.aug,intercept=FALSE,normalize=FALSE,lambda=0)
    iols <- length(senet1$penalty)
    lambda1 <- senet1$penalty[-iols]/2

    ## our approach: direct optimization
    senet2 <- elastic_net(x,y,intercept=intercept,normalize=normalize,
                          struct=t(C) %*% C, lambda1=lambda1, lambda2=lambda2)

    res <- list(
      coef.ref = t(scale(predict(senet1,
        type="coefficients",naive=TRUE)$coefficients,FALSE,normx)[-iols,]),
        coef.our = as.matrix(senet2$coefficients))

    res

  }

  ## PROSTATE DATA SET
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate
  p <- ncol(x)

  ## Simple Elastic.net: structuring matrix is the indentity
  C <- diag(rep(1,p))
  out <- get.coef(x,y,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## Structured Elastic.net
  C <- as.matrix(Matrix::bandSparse(p,k=0:1,diagonals=list(rep(1,p),rep(-1,p-1))))
  ## with intercept and normalization
  out <- get.coef(x,y,intercept=TRUE,normalize=TRUE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## without intercept, with normalization
  out <- get.coef(x,y,intercept=FALSE,normalize=TRUE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## without intercept, without normalization
  out <- get.coef(x,y,intercept=FALSE,normalize=FALSE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## with intercept, without normalization
  out <- get.coef(x,y,intercept=TRUE,normalize=FALSE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
  set.seed(seed)

  beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
  n <- 100
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  ## Simple Elastic.net: structuring matrix is the indentity
  C <- diag(rep(1,p))
  out <- get.coef(x,y,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## Structured Elastic.net
  C <- as.matrix(bandSparse(p,k=0:1,diagonals=list(rep(1,p),rep(-1,p-1))))
  ## with intercept and normalization
  out <- get.coef(x,y,intercept=TRUE,normalize=TRUE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## without intercept, with normalization
  out <- get.coef(x,y,intercept=FALSE,normalize=TRUE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## without intercept, without normalization
  out <- get.coef(x,y,intercept=FALSE,normalize=FALSE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

  ## with intercept, without normalization
  out <- get.coef(x,y,intercept=TRUE,normalize=FALSE,C=C)
  expect_equal(out$coef.our, out$coef.ref, check.attributes = FALSE, tolerance = tol)

})

