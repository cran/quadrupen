context("Consistency of the Elastic-net solution path")

test_that("Consistency between 'quadrupen' and 'elasticnet' packages", {

  require(elasticnet)

  get.coef <- function(x,y,intercept,naive=FALSE) {
    lambda2 <- runif(1,0,10)
    enet.larsen <- enet(x,y,lambda=lambda2,intercept=intercept)
    iols <- length(enet.larsen$penalty)
    lambda1 <- enet.larsen$penalty[-iols]/2
    enet.quadru <- elastic.net(x,y,intercept=intercept,lambda1=lambda1, lambda2=lambda2, naive=naive)
    return(list(coef.quad=as.matrix(enet.quadru@coefficients),
                coef.enet=predict(enet.larsen, type="coefficients",naive=naive)$coefficients[-iols,]))
  }

  ## PROSTATE DATA SET
  prostate <- read.table("http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/prostate.data")
  x <- as.matrix(prostate[,1:8])
  y <- prostate[,9]

  ## Run the tests...
  with.intercept <-get.coef(x,y,intercept=TRUE,naive=TRUE)
  expect_that(with.intercept$coef.quad,
              is_equivalent_to(with.intercept$coef.enet))

  without.intercept <-get.coef(x,y,intercept=FALSE,naive=TRUE)
  expect_that(without.intercept$coef.quad,
              is_equivalent_to(without.intercept$coef.enet))

  with.intercept <-get.coef(x,y,intercept=TRUE,naive=FALSE)
  expect_that(with.intercept$coef.quad,
              is_equivalent_to(with.intercept$coef.enet))

  without.intercept <-get.coef(x,y,intercept=FALSE,naive=FALSE)
  expect_that(without.intercept$coef.quad,
              is_equivalent_to(without.intercept$coef.enet))

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
  with.intercept <-get.coef(x,y,intercept=TRUE)
  expect_that(with.intercept$coef.quad,
              is_equivalent_to(with.intercept$coef.enet))

  without.intercept <-get.coef(x,y,intercept=FALSE)
  expect_that(without.intercept$coef.quad,
              is_equivalent_to(without.intercept$coef.enet))

})

