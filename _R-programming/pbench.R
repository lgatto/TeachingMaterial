library("parallel")
library("microbenchmark")

ll <- replicate(8, matrix(rnorm(1e6),1000), simplify=FALSE)
f <- function(x) mean(solve(x), trim=0.7)

pbench <- microbenchmark(
    res <- lapply(ll, f),
    resmc <- mclapply(ll, f, mc.cores = 16L),
    times = 10)

stopifnot(identical(res, resmc))

save(pbench, file = "pbench.rda")
