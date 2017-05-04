
x <- runif(100)
system.time(sqrt(x))
system.time(x^0.5)


x <- runif(1e5)
system.time(sqrt(x))
system.time(x^0.5)


summary(replicate(10, system.time(x^0.5)[["elapsed"]]))


x <- runif(100)
library(microbenchmark)

microbenchmark(sqrt(x),
               x ^ 0.5)


m <- matrix(rnorm(1e6), ncol = 10)

Rprof("rprof")
res <- apply(m, 1, mean, trim=.3)
Rprof(NULL)
summaryRprof("rprof")


f <- function() {
  pause(0.1)
  g()
  h()
}

g <- function() {
  pause(0.1)
  h()
}

h <- function() {
  pause(0.1)
}



library("profvis")
source("lineprof-example.R")
profvis(f())


profvis({
    m <- matrix(rnorm(1e6), ncol = 10)
    res <- apply(m, 1, mean, trim=.3)
    sum(res)
})


x <- runif(100)
all.equal(sqrt(x), x ^ 0.5)


library("sequences")
gccount
gccountr <- function(x) table(strsplit(x, "")[[1]])
gccountr2 <- function(x) tabulate(factor(strsplit(x, "")[[1]]))


s <- paste(sample(c("A", "C", "G", "T"),
                  100, replace = TRUE),
           collapse = "")

gccount(s)
gccountr(s)
gccountr2(s)


library("microbenchmark")
microbenchmark(gccount(s),
                     gccountr(s),
                     gccountr2(s),
                     times = 1e4,
					 unit = "eps")


library("ggplot2")
mb <- microbenchmark(gccount(s),
                     gccountr(s),
                     gccountr2(s))
print(mb)
microbenchmark:::autoplot.microbenchmark(mb)


make_id2GO <- function(n = 1e3) { ## could be 1e4 - 1e5
    gn <- sprintf(paste0("ENSG%0", 10, "d"), sample(1e6, n))
    goid <- function(n = 10) sprintf(paste0("GO:%0", 10, "d"), sample(1e6, n))
    structure(replicate(n, goid(sample(50, 1))),
              names = gn)
}
id2GO <- make_id2GO()


length(id2GO)
str(head(id2GO))
str(unlist(id2GO))


library(microbenchmark)
microbenchmark(unlist(l),
               unlist(l, use.names = FALSE),
               times = 10)


f1 <- function(n) {
  a <- NULL 
  for (i in 1:n) a <- c(a, sqrt(i))
  a
}

f2 <- function(n) {
  a <- numeric(n)
  for (i in 1:n) a[i] <- sqrt(i)
  a
}


microbenchmark(f1(1e3), f2(1e3))
microbenchmark(f1(1e4), f2(1e4))


e <- new.env()
e$x <- 1
f <- function(myenv) myenv$x <- 2
f(e)
e$x


f3 <- function(n)
  sapply(seq_len(n), sqrt)

f4 <- function(n) sqrt(n)


n <- 10^(2:5)
t1 <- sapply(n, function(.n) system.time(f1(.n))[["elapsed"]])
t2 <- sapply(n, function(.n) system.time(f2(.n))[["elapsed"]])
t3 <- sapply(n, function(.n) system.time(f3(.n))[["elapsed"]])
t4 <- sapply(n, function(.n) system.time(f4(.n))[["elapsed"]])

elapsed <- data.frame(t1, t2, t3, t4)
rownames(elapsed) <- n

colnames(elapsed) <-
    c("for loop\nwithout init",
      "for loop\nwith init",
      "wrapped in\napply",
      "built-in sqrt\n(vectorised)")
	
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("ggplot2"))

mainvp <- viewport(width = 1,
                   height = 1,
                   x = 0.5, y = 0.5)
subvp <- viewport(width = 6/9,
                  height = 5/9,
                  x = .1,
                  y = .95,
                  just = c("left","top"))
df <- melt(elapsed)
colnames(df) <- c("Implementation", "Elapsed")
df$Iterations <- rep(n, 4)
ymax <- max(elapsed[, -1])
p <- ggplot(data=df, aes(x=Iterations, y=Elapsed, col=Implementation)) +
    geom_line() + geom_point() +
        theme(legend.position="bottom") +
            scale_x_continuous(trans=log10_trans()) +
                coord_trans(x="log2")
q <- p + coord_cartesian(ylim=c(0, (ymax+.05))) +
    theme_gray(8) +
        labs(x = NULL, y = NULL) +
            theme(plot.margin = unit(rep(0.3, 4), "lines")) +
                theme(legend.position="none")
print(p, vp = mainvp)
print(q, vp = subvp)


lapply2 <- function(x, f, ...) {
  out <- vector("list", length(x))
  for (i in seq_along(x)) {
    out[[i]] <- f(x[[i]], ...)
  }
  out
}

lapply2_c <- compiler::cmpfun(lapply2)

x <- list(1:10, letters, c(FALSE, TRUE), NULL)


microbenchmark(
  lapply2(x, is.null),
  lapply2_c(x, is.null),
  lapply(x, is.null))


library("pryr")
library("profvis")


x <- 1:1e5
object.size(x)
print(object.size(x), units = "Kb")
object_size(x)


ll <- list(x, x, x)
print(object.size(ll), units = "Kb")
object_size(ll)


x <- 1:1e6
y <- list(1:1e6, 1:1e6, 1:1e6)
object_size(x)
object_size(y)


e <- new.env()
object.size(e)
object_size(e)
e$a <- 1:1e6
object.size(e)
object_size(e)


mem_used()


mem_change(v <- 1:1e6)
mem_change(rm(v))


rm(list = ls())
mem_change(x <- 1:1e6)
mem_change(y <- x)
mem_change(rm(x))
mem_change(rm(y))

