
```r
library("parallel")

n <- 10000
f1 <- function(n) {
    l <- list()
    for (i in seq_len(n)) l[[i]] <- seq(i)
    return(l)
}

f2 <- function(n) {
    l <- vector("list", length = n)
    for (i in seq_len(n)) l[[i]] <- seq(i)
    return(l)
}

f3 <- function(n) lapply(seq_len(n), seq)


f4 <- function(n, nc = 2L) mclapply(seq_len(n), seq, mc.cores = nc)

library("rbenchmark")
benchmark(f1(n), f2(n), f3(n), f4(n), columns = c("test", "replications", "elapsed", 
    "relative"), replications = 10)
```

```
##    test replications elapsed relative
## 1 f1(n)           10  18.023    3.034
## 2 f2(n)           10   5.941    1.000
## 3 f3(n)           10   9.368    1.577
## 4 f4(n)           10  28.531    4.802
```

