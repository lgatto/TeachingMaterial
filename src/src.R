## Example script file, to generate a little plot.
## Demonstrate trigonometric functions.
## Sept 2007

## @knitr trig
x <- seq(from=0, to=2*pi, length=100)
y <- sin(x)
z <- cos(2*x)
z ## will not appear when source'd
print(y[1:10]) ## should use print()
plot(x, y, type='l')
lines(x, z, type='l', col='red')

## @knitr simplenorm
#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

stopifnot(length(args)==3)
args = as.numeric(args)
n = args[1]; mean = args[2]; sd = args[3]
rnorm(n, mean, sd)


## @knitr mfrow
pdf(file = 'mfrow_eg.pdf',
    width = 4, height = 6)
par(mfrow = c(3,2))
par(mar = c(3.5, 3.5, 1.5, 0.5),
    mgp = c(2.5, 1, 0))
x <- seq(from = 0, to = 2*pi,
         len = 100)
plot(x, sin(x), main = "sin (x)",
     type = 'l')
plot(x, sin(2*x), main = "sin (2x)",
     type = 'l')
plot(x, sin(3*x), main = "sin (3x)",
     type = 'l')
plot(x, cos(x), main = "cos (x)",
     type = 'l')
plot(x, cos(2*x), main = "cos (2x)",
     type = 'l')
plot(x, cos(3*x), main = "cos (3x)",
     type = 'l')
dev.off()

## @knitr debug1
start <- function() { go( sqrt(10)) }
go <- function(x) { inner(x, '-13')}
inner <- function(a, b) {
  c <- sqrt(b)
  a * log(c)
}
start() ## error
traceback() ## postmortem debugging
