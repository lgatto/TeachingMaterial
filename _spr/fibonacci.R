### R code from vignette source 'fibonacci.Rnw'

###################################################
### code chunk number 1: external-src
###################################################
library(knitr)
options(width = 60)
opts_chunk$set(prompt = TRUE,
               comment = '',
               fig.align = 'center')



###################################################
### code chunk number 2: fibo1
###################################################
fibonacci <- function(n) {
  .fib <- function(n) {
    if (n == 0)
      return(0)
    if (n < 3)
      return(1)
    ans <- c(1, 1, rep(0, n-2))
    for (i in 3:n)
      ans[i] <- ans[i-1] + ans[i-2]
    return(ans[n])
  }
  if (length(n) == 1)
    return(.fib(n))
  sapply(n, .fib)
}


###################################################
### code chunk number 3: fibodirect
###################################################
fibdirect <- function(n) {
  stopifnot(n >= 0)
  phi <- (1+sqrt(5))/2 
  round((phi^n - (1 - phi)^n)/sqrt(5))
}


###################################################
### code chunk number 4: checkfib
###################################################
x <- 1:100
all.equal(fibonacci(x), fibdirect(x))


###################################################
### code chunk number 5: fibbench
###################################################
library("rbenchmark")
benchmark(fibonacci(x), fibdirect(x),
          columns=c("test", "replications", 
            "elapsed", "relative"),
          order = "relative", replications = 100)


###################################################
### code chunk number 6: gratio
###################################################
all.equal((1+sqrt(5))/2, fibonacci(10)/fibonacci(9))
all.equal((1+sqrt(5))/2, fibonacci(20)/fibonacci(19))
all.equal((1+sqrt(5))/2, fibonacci(30)/fibonacci(29))


###################################################
### code chunk number 7: fibplot
###################################################
phi <- (1+sqrt(5))/2
plot(fibdirect(2:20)/fibdirect(1:19), type = "b")
abline(h = phi, lty = "dotted", col = "red")
plot(fibdirect(19:27)/fibdirect(18:26), type = "b")
abline(h = phi, lty = "dotted", col = "red")


