
f <- function(x) {
    y <- x + 1
    return(x * y)
}


body(f)
args(f)
environment(f)

body(f) <- quote({
    y <- x * y
    return(x + y)
})


f <- function(x) x + y

f(1)

environment(f)
y <- 2
f(1)



e <- new.env()
environment(f) <- e

f(1)
e$y <- 10
f(1)


codetools::findGlobals(f)


f <- function() {
    x <- 1
    y <- 2
    c(x, y)
}
f()


x <- 2
g <- function(){
    y <- 1
    c(x, y)
}
g()


x <- 1
h <- function() {
    y <- 2
    i <- function() {
        z <- 3
        c(x, y, z)
    }
    i()
}
h()


j <- function(x) {
    y <- 2
    function(){
        c(x, y)
    }
}
k <- j(1)
k()


j <- function() {
    if (!exists("a")) {
        a <- 1
    } else {
        a <- a + 1
    }
    print(a)
}
j() ## First call
j() ## Second call


f <- function() x
x <- 1
f()
x <- 2
f()


f <- function(x) {
    f <- function(x) {
        f <- function(x) {
            x^2
        }
        f(x) + 1
    }
    f(x) * 2
}
f(10)


args <- list(x = 1:10, trim = 0.3)
do.call(mean, args)


f <- function(x = 1, y = 2) x * y
f <- function(x = 1, y = x + 2) x * y


f <- function(x = 1, y) {
	c(missing(x), missing(y))
}
f()
f(x = 1)


plot2 <- function(...) {
    message("Verbose plotting...")
    plot(...)
}

f <- function(...) list(...)


f1 <- function() 1
f2 <- function() return(1)
f3 <- function() return(invisible(1))


f1 <- function(x) {
    on.exit(print("!"))
    x + 1
}

f2 <- function(x) {
    on.exit(print("!"))
    stop("Error")
}


f3 <- function() {
    on.exit(print("1"))
    on.exit(print("2"))
    invisible(TRUE)
}


f4 <- function() {
    on.exit(print("1"))
    on.exit(print("2"), add = TRUE)
    invisible(TRUE)
}


function(x) x + y
body(function(x) x + y)
args(function(x) x + y)
environment(function(x) x + y)


make.power <- function(n)
    function(x) x^n


cube <- make.power(3)
square <- make.power(2)
cube(2)
square(2)
environment(cube)
environment(square)


new_counter <- function() {
    i <- 0
    function() {
        i <<- i + 1
        i
    }
}

count1 <- new_counter()
count2 <- new_counter()

count1()
count1()
count2()

environment(count1)
environment(count2)
environment(count1)$i
environment(count2)$i


colramp <- colorRampPalette(c("blue", "yellow"))
colramp(5)
plot(1:10, col = colramp(10), pch = 19, cex = 2,
     main = "colramp(10)")


L <- replicate(3, matrix(rnorm(9), 3), simplify = FALSE)
Reduce("+", L)
try(sum(L))


Reduce("+", list(1, 2, 3), init = 10)
Reduce("+", list(1, 2, 3), accumulate = TRUE)
Reduce("+", list(1, 2, 3), right = TRUE, accumulate = TRUE)


even <- function(x) x %% 2 == 0
(y <- sample(100, 10))
Filter(even, y)
Filter(Negate(even), y)


Map(even, 1:3)


Find(even, 10:15)
Find(even, 10:15, right = TRUE)
Position(Negate(even), 10:15)
Position(Negate(even), 10:15, right = TRUE)


(x <- 1:5)
(y <- 5:1)
x + y


(x <- 1:6)
(y <- 1:2)
x+y


diff1 <- function(e) {
  n <- length(e)
  interval <- rep(0, n - 1)
  for (i in 1:(n - 1))
      interval[i] <- e[i + 1] - e[i]
  interval
}
e <- c(2, 5, 10.2, 12, 19)
diff1(e)


diff2 <- function(e) {
  n <- length(e)
  e[-1] - e[-n]
}
e <- c(2, 5, 10.2, 12, 19)
diff2(e)


v <- rnorm(1000) ## or a list
res <- numeric(length(v))

for (i in 1:length(v))
  res[i] <- f(v[i])

res <- sapply(v, f)

## if f is vectorised
f(v)


M <- matrix(rnorm(100), 10)
apply(M, 1, function(Mrow) 'do something with Mrow')
apply(M, 2, function(Mcol) 'do something with Mcol')


f <- function(x, a = 1) sin(x^2)/ (a + abs(x))
x <- seq(-7, 7, 0.02 )
x0 <- seq(-2, 2, 0.02)
y0 <- f(x0)
y0[y0 < 0] <- 0
plot(x, f(x), type = "l", main = expression(f(x) ==  frac(sin(x^2),(a + abs(x)))))
grid()
abline(v = c(-2, 2), lty = "dotted")
polygon(x0, y0, col = "#00000010")


f <- function(x, a = 1) sin(x^2)/ (a + abs(x))
integrate(f, lower = -2, upper = 2)


lo <- c(-2, 0)
hi <- c(0, 2)
integrate(f, lower = lo, upper = hi)


mapply(function(lo, hi) integrate(f, lo, hi)$value,
       lo, hi)


Integrate <- Vectorize(
  function(fn, lower, upper)
  integrate(fn, lower, upper)$value,
  vectorize.args=c("lower", "upper")
  )
Integrate(f, lower=lo, upper=hi)

