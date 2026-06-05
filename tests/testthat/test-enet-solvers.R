context("Testing consistency and timings of the different solvers for elastic-net")

testDataEnet <- readRDS("dataTest-Enet.rds")

test_that("Cholesky, conjugate gradient, FISTA, PGD ", {

  tol <- 1e-2
  
  get.coef <- function(x,y) {
    lambda1 <- .25

    enet.ref <- elastic_net(x,y,lambda1=lambda1, control=list(timer=TRUE))
    enet.cg  <- elastic_net(x,y,lambda1=lambda1, control=list(timer=TRUE,usechol=FALSE))
    enet.fista <- elastic_net(x,y,lambda1=lambda1, control=list(timer=TRUE,method = "fista"))
    enet.pgd <- elastic_net(x,y,lambda1=lambda1, control=list(timer=TRUE,method = "pgd"))

    cat("\n\tTimings with warm-restart along the path")
    cat("\n\t\tQuadratic from stratch (Cholesky): ",enet.ref$optim_monitoring$timer)
    cat("\n\t\tQuadratic from stratch (Conjugate Gradient): ", enet.cg$optim_monitoring$timer)
    cat("\n\t\tFISTA from stratch: ", enet.fista$optim_monitoring$timer)
    cat("\n\t\tPGD + ANDERSON from stratch: ", enet.pgd$optim_monitoring$timer)

    return(list(
      coef.ref=as.matrix(enet.ref$coefficients),
      coef.cg =as.matrix(enet.cg$coefficients),
      coef.fista =as.matrix(enet.fista$coefficients),
      coef.pgd =as.matrix(enet.pgd$coefficients)
      )
    )
      
  }

  ## PROSTATE DATA SET
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate
  
  ## Run the tests...
  cat("\n  * tiny-size problem...")
  out <- get.coef(x,y)
  cat("\n")
  expect_equal(out$coef.cg   ,out$coef.ref    ,tolerance=tol)
  expect_equal(out$coef.fista,out$coef.ref    ,tolerance=tol)
  expect_equal(out$coef.pgd  ,out$coef.ref    ,tolerance=tol)

  ## RANDOM DATA
  x <- testDataEnet$x_sim
  y <- testDataEnet$y_sim
  
  ## Run the tests...
  cat("\n  * small-size problem, with correlation...")
  out <- get.coef(x,y)
  cat("\n")
  expect_equal(out$coef.cg   ,out$coef.ref    ,tolerance=tol)
  expect_equal(out$coef.fista,out$coef.ref    ,tolerance=tol)
  expect_equal(out$coef.pgd  ,out$coef.ref    ,tolerance=tol)
  
})

