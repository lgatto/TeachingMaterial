find.high <- function(x, t) {
  ## Return samples in x bigger than t.
  ## (Better to use x[x>t] in real-life!)
  max.length <- 100  ## should be upper limit...
  results <- rep(0, max.length)
  counter <- 0
  for (i in x) {
    if (i > t) {
      counter <- counter + 1
      if (counter > max.length) {
        browser()
      } else {
        results[counter] <- i
      }
    }
  }
  results[1:counter]
}
x <- rnorm(100)
find.high(x, 0.7)

x <- rnorm(1000)
(1- pnorm(0.7)) * length(x) ## expected.
find.high(x, 0.7)

