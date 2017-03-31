
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


f <- function() 1
f()
body(f)
body(f) <- 2
f()
formals(f) <- pairlist(x = 1)
body(f) <- quote(x + 1)
f()
f(10)


setMethod(...)


ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
mod <- lm(weight ~ group)
mod
mod$foo <- 1
mod
summary(mod)
names(mod)


class(mod) <- c("mod", "foo")


## function
f <- function(x) NULL

## S3 method
s3 <- function(x) UseMethod("s3")
s3.integer <- f

## S4 method
A <- setClass("A", slots = c(a = "list"))
setGeneric("s4", function(x) standardGeneric("s4"))
setMethod("s4", "A", f)
a <- A()

## S4 Reference Class
B <- setRefClass("B", methods = list(rc = f))
b <- B$new()


microbenchmark(
     fun = f(),
     S3 = s3(1L),
     S4 = s4(a),
     RC = b$rc()
)


a <- 1
f <- function() {
    g <- function() {
        print(a) ## from global workspace
        assign("a", 2, envir = parent.frame())
        print(a) ## f's environment
        a <- 3
        print(a) ## g's environment
    }
    g()
}
f()


f <- function(x, y) {
    (x + y) ^ 2
}


microbenchmark(
    "[32, 11]"      = mtcars[32, 11],
    "$carb[32]"     = mtcars$carb[32],
    "[[c(11, 32)]]" = mtcars[[c(11, 32)]],
    "[[11]][32]"    = mtcars[[11]][32],
    ".subset"       = .subset(mtcars, 11)[32],
    ".subset2"      = .subset2(mtcars, 11)[32]
)


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


x <- 1:10
c(address(x), refs(x))
x[5] <- 0L
c(address(x), refs(x))


y <- x
c(address(x), refs(x))
c(address(y), refs(y))

x[5] <- 1L
c(address(x), refs(x))
c(address(y), refs(y))


x <- 1:5
y <- x
c(address(x), refs(x))
rm(y)
c(address(x), refs(x)) ## should be 1

x <- 1:5
y <- x
z <- x
c(address(x), refs(x)) ## should be 3


x <- 1:10
tracemem(x)
x[5] <- 0L

y <- x
x[5] <- 10L

address(x)
address(y)


f1 <- function(n) {
  a <- NULL
  for (i in 1:n) {
      a <- c(a, sqrt(i))
      print(address(a))
  }
  invisible(a)
}

f2 <- function(n) {
  a <- numeric(n)
  for (i in 1:n) {
      a[i] <- sqrt(i)
      print(address(a))
  }
  invisible(a)
}

