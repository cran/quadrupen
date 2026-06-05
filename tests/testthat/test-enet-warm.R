context("Consistency and timings of warm restart for the elastic_net")

testDataEnet <- readRDS("dataTest-Enet.rds")

test_that("Warm restart works for Elastic-Net", {

  require(quadrupen)

  get.coef <- function(x,y) {
    lambda1 <- .25
    
    enet.ref <- elastic_net(x,y,lambda1=lambda1, control=list(timer=TRUE))

    enet.ref.bot <- elastic_net(x,y,lambda1=lambda1*2)
    enet.ref.up  <- elastic_net(x,y,lambda1=lambda1/2)
    beta0_bt <- as.numeric(enet.ref.bot$coefficients)
    beta0_up <- as.numeric(enet.ref.up$coefficients)
    enet.bot <- elastic_net(x,y,lambda1=lambda1, beta0 = beta0_bt, control=list(timer=TRUE))
    enet.up  <- elastic_net(x,y,lambda1=lambda1, beta0 = beta0_up, control=list(timer=TRUE))

    cat("\n\tTimings with warm-restart along the path")
    cat("\n\t\tfrom stratch: ",enet.ref$optim_monitoring$timer)
    cat("\n\t\tstarting from sparser solution: ",enet.bot$optim_monitoring$timer)
    cat("\n\t\tstarting from more dense solution: ",enet.up$optim_monitoring$timer)
    cat("\n\n")
    
    return(list(
      coef.ref=as.matrix(enet.ref$coefficients),
      coef.bot=as.matrix(enet.bot$coefficients),
      coef.up =as.matrix(enet.up$coefficients)))
  }

  ## PROSTATE DATA SET
  x <- testDataEnet$x_prostate
  y <- testDataEnet$y_prostate
  
  ## Run the tests...
  cat("\n  * tiny-size problem...")
  out <- get.coef(x,y)
  expect_equal(out$coef.bot, out$coef.ref, check.attributs = FALSE)
  expect_equal(out$coef.up , out$coef.ref, check.attributs = FALSE)

})

