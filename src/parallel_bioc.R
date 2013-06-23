library(BiocParallel)

ll <- replicate(8, matrix(rnorm(1e6),1000), simplify=FALSE)
f <- function(x) mean(solve(x), trim=0.7)

(p <- SnowParam(4L))
bplapply(ll, f, BPPARAM = p)

