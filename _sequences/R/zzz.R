.onAttach <- function(libname, pkgname) {
  msg <- paste("This is package 'sequences'\n")
  packageStartupMessage(msg)
  ## addVigs2WinMenu("sequences") ## requires Biobase
}
