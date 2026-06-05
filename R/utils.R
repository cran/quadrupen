chain_graph <- function(p) {
  G <- list(conn = list(), weight = list())
  G$conn[1:(p - 1)] <- as.integer(2:p) 
  G$conn[p] <- as.integer(1)
  G$weight[1:p] <- rep(1.0, p)
  class(G) <- "FusionGraph"
  G
}

conn_from_adj <- function(adjacency_matrix) {
  G <- list(conn = list(), weight = list())
  pairs <- which(upper.tri(adjacency_matrix), arr.ind = TRUE)
  for (idx in 1:nrow(pairs)) {
    i = unname(pairs[idx,][1])
    j = unname(pairs[idx,][2])
    if (adjacency_matrix[i, j] == 1) {
      G$conn[[i]] = j
      G$conn[[j]] = i
      G$weight[[i]] = G$weight[[j]] <- as.numeric(1)
    }
  }
  class(G) <- "FusionGraph"
  G
}

optim_enet_default <- function(d)
  list(verbose     = 0, # default control options
       timer       = FALSE,
       maxiter     = 50,
       method      = "quadra",
       threshold   = 1e-6,
       monitor     = 0,
       factmat     = TRUE
  )

optim_breg_default <- function(d)
  list(verbose     = 0, # default control options
       timer       = FALSE,
       maxiter     = 10,
       method      = "quadra",
       threshold   = 1e-4,
       monitor     = 0,
       factmat     = FALSE
  )

optim_fused_default <- function(d)
  list(verbose       = 0, # default control options
       timer         = FALSE,
       maxiterout    = 100,
       maxiterin     = 10000,
       maxactivation = 10, 
       accuracy      = 1e-6,
       monitor       = 0,
       fusioncheck   = "all" ## c("all","active","none", "naive", "huber")
  )

optim_grp_default <- function(d)
  list(verbose     = 0, # default control options
       timer       = FALSE,
       maxiter     = 50,
       method      = "fista",
       threshold   = 1e-3,
       monitor     = 0,
       factmat     = FALSE
  )

status_to_message <- function(status) {
  message <- switch(as.character(status),
                    "0"  = "converged",
                    "1"  = "max # of iterate reached",
                    "2"  = "max # of feature reached",
                    "3"  = "system has become singular",
                    "Return status not recognized"
  )
  message
}


#' @importFrom rlang .data
.plot_regpath <- function(xv, coef, standardize, log_scale, labs) {

  nzeros <- attr(coef, "nzeros")  
  
  dplot <- data.frame(xvar = xv, coef = coef) |> 
    tidyr::pivot_longer(cols = -.data$xvar, names_to = "var", values_to = "coef")
  
  if (is.null(labs)) {
    dplot$labs <- factor(rep(nzeros, length(xv)))
  } else {
    if (sum(is.na(labs[nzeros])) > 0 ) {
      labs <- NULL
      warning("The number of label is wrong, ignoring them.")
      dplot$labs <- factor(rep(nzeros, length(xv)))
    } else {
      dplot$labs <- factor(rep(labs[nzeros], length(xv)))
    }
  }
  
  d <- ggplot(dplot) + aes(x = .data$xvar, y = coef, color = labs, group = .data$var) + 
    geom_line() +  geom_hline(yintercept = 0, alpha = 0.5, linetype = "dotted") +
    ylab(ifelse(standardize, "standardized coefficients","coefficients"))

  if (attr(xv, "type") == "lambda") {
    d <- d + xlab(ifelse(log_scale,expression(log[10](lambda)),expression(lambda)))
    if (log_scale)
      d <- d + scale_x_log10() + annotation_logticks(sides = "b")
  } else if (attr(xv, "type") == "fraction") {
    d <- d + xlab(expression(paste("|",beta[lambda],"|",{}[1]/max[lambda],"|",beta[lambda],"|",{}[1],sep = "")))
  } else {
    d <- d + xlab("Degrees of freedom")
  }
  d
}
