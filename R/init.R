.onAttach <- function(...) {
  welcome <- paste(""                                                          ,
                   "---------------------------------------------------",
                   ""                                                    ,
                   "  'quadrupen' package version 0.2-0                " ,
                   ""                                                    ,
                   " May be buggy and default parameters inappropriate" ,
                   " Any user's feedback is thus very welcome!",
                   "---------------------------------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

