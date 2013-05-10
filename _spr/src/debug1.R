start <- function() { go( sqrt(10)) }
go <- function(x) { inner(x, '-13')}
inner <- function(a, b) {
  c <- sqrt(b)
  a * log(c)
}
start()
traceback()

