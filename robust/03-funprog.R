
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


x <- 1
i <- function() {
    z <- 3
    c(x, y, z)
}
h <- function() {
    y <- 2
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


v <- rnorm(1000) ## or a list
res <- numeric(length(v))

for (i in 1:length(v))
  res[i] <- f(v[i])

res <- sapply(v, f)

## if f is vectorised
f(v)


sqrtabs <- function(x) {
    v <- abs(x)
    sapply(1:length(v), function(i) sqrt(v[i]))
}


M <- matrix(rnorm(100), 10)
apply(M, 1, function(Mrow) 'do something with Mrow')
apply(M, 2, function(Mcol) 'do something with Mcol')


df1 <- data.frame(x = 1:3, y = LETTERS[1:3])
sapply(df1, class)
df2 <- data.frame(x = 1:3, y = Sys.time() + 1:3)
sapply(df2, class)


lapply(df1, class)
lapply(df2, class)


vapply(df1, class, "1")
vapply(df2, class, "1")

