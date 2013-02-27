.onAttach <- function(...) {
  welcome <- paste(""                                              ,
                   "----------------------------------------------",
                   "  'quadrupen' package version 0.2-2           ",
                   ""                                              ,
                   " Still under development... feedback welcome  ",
                   "----------------------------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

