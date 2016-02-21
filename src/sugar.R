source("src/pi.R")

library("Rcpp")
sourceCpp("src/pi.cpp")

N <- 1e6
set.seed(42)
resR <- piR(N)

set.seed(42)
resCpp <- piSugar(N)
stopifnot(identical(resR, resCpp))


library(rbenchmark)
res <- benchmark(piR(N),
                 piSugar(N),
                 order="relative")
print(res[,1:4])
