.onAttach <- function(...) {
  welcome <- paste("'quadrupen' package version 1.0-0", sep = "\n")
  packageStartupMessage(welcome)
}

