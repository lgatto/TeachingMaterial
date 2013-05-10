## Example for demonstrating debug()
## The "if ... browser ..." line can be deleted and added later.
## source('stddev.R')

std.dev <- function(x) {
  ## Return std dev of X.
  n <- length(x)
  if (n <= 1) {
    browser()
  }
  xbar <- sum(x)/n
  browser()
  diff <- x - xbar
  sum.sq <- sum( diff^2 )
  var <- sum.sq / (n-1)
  ans <- sqrt(var)
  return(ans)
}



x <- seq(1:9)
std.dev(x)
             
