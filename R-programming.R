
## @knitr env, include=FALSE, echo=FALSE, cache=FALSE
library("knitr")
opts_chunk$set(fig.align = 'center', 
               fig.show = 'hold', 
               tidy = FALSE,
               par = TRUE,
               prompt = TRUE,
               eval = TRUE,
               stop_on_error = 1L,
               comment = NA)
options(replace.assign = TRUE, 
        width = 55)
set.seed(1)


## @knitr fun, tidy = FALSE
myFun <- function(i, j = 1) {
    mn <- min(i, j)
    mx <- max(i, j)    
    k <- rnorm(ceiling(i * j))
    return(k[k > mn/mx])
}
myFun(1.75, 4.45)
myFun(1.75) ## j = 1 by default


## @knitr scoping
x <- 1
f1 <- function(x) {
    x <- x + 10
    x
}

f1(x)
x ## unchanged


## @knitr scoping2
f2 <- function() {
    x <- x + 10
    x
}

f2()
x ## still unchanged


## @knitr fun2fun
make.power <- function(n)
    function(x) x^n
square <- make.power(2)
cube <- make.power(3)


## @knitr fun2funexplore
square
get("n", environment(square))
square(2)
cube(2)


## @knitr fun2fun2
(rbramp <- colorRampPalette(c("red", "blue")))
rbramp(3)
rbramp(7)


## @knitr plt
plt <- function(n, ...)
    plot(1:n, ...)


## @knitr pltex, dev='pdf', fig.width = 8, fig.height = 4
par(mfrow = c(1, 2))
plt(5, pch = 19, type = "b")
plt(10, col = rbramp(10), pch = 15)


## @knitr args
args(cat)
args(rm)


## @knitr lapply
lapply(1:2, rnorm)
lapply(1:2, rnorm, 10, 2)


## @knitr sapply
library(fortunes)
lapply(sample(315, 1), fortune)
sapply(sample(315, 1), fortune)


## @knitr apply
set.seed(10)
m <- matrix(rnorm(10), ncol = 2)
apply(m, 1, myFun)
apply(m, 1, myFun)
apply(m, 1, max) ## Biobase::rowMax
apply(m, 2, min) ## Biobse::rowMin


## @knitr mapply
mapply(rep, 1:4, 4:1)


## @knitr tapply
dfr <- data.frame(f1 = sample(LETTERS[1:2], 10, replace = TRUE),
                  f2 = sample(LETTERS[3:4], 10, replace = TRUE),
                  x = rnorm(10))
tapply(dfr$x, dfr$f1, mean)
tapply(dfr$x, dfr$f2, mean)
tapply(dfr$x, list(dfr$f1, dfr$f2), mean)


## @knitr anon
m
apply(m, 1, function(x) ifelse(mean(x) > 0, mean(x), max(x)))


## @knitr N, echo = FALSE
N <- 1e4


## @knitr ll, cache = TRUE
ll <- lapply(sample(N), rnorm)
f <- function(x) mean(x) * length(x)


## @knitr time1, cache = TRUE
res1 <- c()
system.time({
    for (i in 1:length(ll))
        res1[i] <- f(ll[[i]])
})


## @knitr time2, cache = TRUE
res2 <- numeric(length(ll))
system.time({
    for (i in 1:length(ll))
        res2[i] <- f(ll[[i]])
})


## @knitr time3, cache = TRUE
system.time(res3 <- sapply(ll, f))


## @knitr replicatetime, cache = TRUE
summary(replicate(50, system.time(res3 <- sapply(ll, f))["elapsed"]))


## @knitr benchmarkfun
sol2 <- function(x) {
    n <- length(x)
    ans <- numeric(n)
    for (i in 1:n) {
        ans[i] <- f(x[[i]])
    }
    ans
}
sol3 <- function(x)
    sapply(x, f)    


## @knitr benmark, cache = TRUE
library("microbenchmark")
microbenchmark(sol2(ll), sol3(ll), times = 200)


## @knitr profiling, cache = TRUE
Rprof()
tmp <- replicate(10, sol3(ll))
Rprof(NULL)


## @knitr opts0, echo=FALSE
oldwidth <- options()$width
options(width = 100)


## @knitr smryprof, size = 'small'
summaryRprof()


## @knitr opts1, echo=FALSE
options(width = oldwidth)


## @knitr id1
identical(res1, res2)


## @knitr id3
identical3 <- 
    function(x,y,z) identical(x,y) && identical (y,z)
identical3(res1, res2, res3)


## @knitr sqrtx
x <- sqrt(2)
x * x == 2
identical(x*x, 2)


## @knitr alleqsqrt
all.equal(x * x, 2)


## @knitr stopifnot
stopifnot(x * x == 2)
stopifnot(all.equal(x * x, 2))


## @knitr pvec
library("parallel")
detectCores()
mclapply(1:3, function(x) Sys.getpid(), mc.cores = 3)
mclapply(1:3, function(x) Sys.getpid(), mc.cores = 2)


## @knitr solpar, cache = TRUE
solmc <- function(x)
    mclapply(x, f)
solpar <- function(x, cl)
    parLapply(cl, x, f)
sol3 <- function(x)
    lapply(x, f)
cl <- makeCluster(4)
stopifnot(identical3(sol3(ll), solmc(ll), solpar(ll, cl)))
stopCluster(cl)


## @knitr pbench, echo = FALSE
cat(scan('pbench.R', what = "", strip.white = FALSE, sep = "\n"), sep = "\n")


## @knitr printpbench, echo=FALSE
load("pbench.rda")
microbenchmark:::print.microbenchmark(pbench)


## @knitr oops, echo=FALSE
e <- function(i) {
  x <- 1:4
  if (i < 5) x[1:2]
  else x[-1:2]
}
f <- function() sapply(1:10, e)
g <- function() f()


## @knitr error, eval=FALSE, prompt = FALSE
## > g()
## Error in x[-1:2] (from #3) : only 0's may be mixed with negative subscripts
## > g
## function() f()


## @knitr traceback, eval=FALSE, prompt = FALSE
## > traceback()
## 5: FUN(1:10[[5L]], ...)
## 4: lapply(X = X, FUN = FUN, ...)
## 3: sapply(1:10, e) at #1
## 2: f() at #1
## 1: g()


## @knitr erroronly, eval=FALSE
## Error in x[-1:2] (from #3) : only 0's may be mixed with negative subscripts


## @knitr showe, eval=FALSE, prompt = FALSE
## e
## function(i) {
##   x <- 1:4
##   if (i < 5) x[1:2]
##   else x[-1:2]
## }
## e(5)
## Error in x[-1:2] (from #3) : only 0's may be mixed with negative subscripts


## @knitr debugmode, eval=FALSE, prompt = FALSE
## > debug(e)
## > e(5)
## debugging in: e(5)
## debug at #1: {
##     x <- 1:4
##     if (i < 5)
##         x[1:2]
##     else x[-1:2]
## }
## Browse[2]>
## debug at #2: x <- 1:4
## Browse[2]>
## debug at #3: if (i < 5) x[1:2] else x[-1:2]
## Browse[2]> ls()
## [1] "i" "x"
## Browse[2]> i
## [1] 5
## Browse[2]> x
## [1] 1 2 3 4
## Browse[2]>
## debug at #3: x[-1:2]
## Browse[2]> x[-1:2]
## Error in x[-1:2] (from #3) : only 0's may be mixed with negative subscripts
## Browse[2]> x[-(1:2)]
## [1] 3 4
## Browse[2]> Q
## > undebug(e)
## > fix(e)


## @knitr sessioninfo, results='asis', echo=FALSE
toLatex(sessionInfo())


