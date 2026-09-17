.onAttach <- function(...) {
  welcome <- paste("'quadrupen' package version 1.1-0", sep = "\n")
  packageStartupMessage(welcome)
}

