##' A simple function for the QuickPackage demo.
##'
##' The function calls the \code{\link{packageDescription}}
##' function to retrieve the code{QuickPackage} description. 
##' @title The QuickPackage description
##' @return An object of class \code{packageDescription}.
##' @author Laurent Gatto
##' @examples
##' qpf()
qpf <-
function() 
  packageDescription("QuickPackage")
