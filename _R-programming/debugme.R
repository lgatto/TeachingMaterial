e <- function(i) {
  x <- 1:4
  if (i < 5) x[1:2]
  else x[-1:2]
}

f <- function() sapply(1:10, e)

g <- function() f()
