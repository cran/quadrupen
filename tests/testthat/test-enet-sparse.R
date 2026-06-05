context("Consistency between sparse/non-sparse encoding")

require(Matrix)
        
test_that("Consistency of quadrupen between sparse/non-sparse encoding of the predictors", {

  ## data generation
  rlm.sparse <- function(n, p, reg="low", mu=3, size=5, prob=0.01) {
    s <- switch(reg,
                low  = floor(0.50 * min(n,p)),
                med  = floor(0.10 * min(n,p)),
                high = floor(0.01 * min(n,p)))

    w <- rep(0,p)
    w[sample(1:p,s)] <- sample(c(-1,1),s,replace=TRUE)*runif(s,1,2)

    X <- matrix(rbinom(n*p, size, prob),n,p)

    Xw <- crossprod(t(X),w)
    sigma <- sum(Xw^2) * 0.01/n
    epsilon <- rnorm(n) * sigma
    y <- mu + Xw + epsilon
    r2 <- 1 - sum(epsilon^2) / sum((y-mean(y))^2)

    return(list(y=y, x=X, r2=r2, w=w, mu=mu))
  }

  n <- 500
  p <- 10000
  lambda2 <- 0.05
  data <- rlm.sparse(n, p, prob=0.01, reg="med")
  s <- sum(data$w != 0)
  y    <- data$y
  x.ns <- as.matrix(data$x)
  x.sp <- Matrix(data$x, sparse=TRUE)
  max.feat <- s * 10

  cat("\n\tProblem with",p, "predictors and",n,"samples,",s,"true nonzeros in beta.star.")
  cat("\n\tThe densely encoded design matrix weights" , round(object.size(x.ns) /(1024^2),2), "Mo.")
  cat("\n\tThe sparsely encoded design matrix weights", round(object.size(x.sp) /(1024^2),2), "Mo.")

  cat("\n\tdense coding...")
  out.enet.ns <- elastic_net(x.ns, y, lambda2=lambda2, maxfeat=max.feat, minratio=1e-3, control=list(timer=TRUE))

  cat(" took", out.enet.ns$optim_monitoring$timer, "seconds to activate",
      colSums(out.enet.ns$coefficients!=0)[length(out.enet.ns$major_tuning)],"variables.")

  cat("\n\tsparse coding...")
  out.enet.sp <- elastic_net(x.sp, y, lambda2=lambda2, maxfeat=max.feat, minratio=1e-3, control=list(timer=TRUE))
  cat(" took", out.enet.sp$optim_monitoring$timer, "seconds to activate",
      colSums(out.enet.sp$coefficients!=0)[length(out.enet.sp$major_tuning)],"variables.\n")
  
  ## remove monitoring for fair comparison!!!
  out.enet.sp$optim_monitoring <- list()
  out.enet.ns$optim_monitoring <- list()
  expect_equal(out.enet.sp$coefficients, out.enet.ns$coefficients, check.attributes = FALSE, tolerance = 1e-4)
  expect_equal(out.enet.sp$interceptTerm, out.enet.ns$interceptTerm, check.attributes = FALSE, tolerance = 1e-4)
})


