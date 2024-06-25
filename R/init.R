.onAttach <- function(...) {
  welcome <- paste(""                                                      ,
                   "------------------------------------------------------",
                   "  'quadrupen' package version 0.2-12                   ",
                   ""                                                      ,
                   " Dev version on https://github.com/jchiquet/quadrupen ",
                   "------------------------------------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

