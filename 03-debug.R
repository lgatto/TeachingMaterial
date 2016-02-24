
message("This is a message for our dear users.")


message("This is a message for our dear users. ",
	paste("Thank you for using our software",
              sw, "version", packageVersion(sw)))


f1 <- function() {
    cat("I AM LOUD AND YOU CAN'T HELP IT.\n")
    ## do stuff
    invisible(TRUE)
}
f1()


f2 <- function() {
    message("Sorry to interup, but...")
    ## do stuff
    invisible(TRUE)
}
f2()
suppressMessages(f2())


f3 <- function(verbose = TRUE) {
    if (verbose)
        message("I am being verbose because you let me.")
    ## do stuff
    invisible(TRUE)
}
f3()
f3(verbose = FALSE)


warning("Do not ignore me. Somthing bad might have happened.")
warning("Do not ignore me. Somthing bad might be happening.", immediate. = TRUE)


f <- function(...)
    warning("Attention, attention, ...!", ...)
f()
f(call. = FALSE)


warnings()
last.warning


option("warn")


stop("This is the end, my friend.")


log(c(2, 1, 0, -1, 2)); print('end') ## warning
xor(c(TRUE, FALSE));  print ('end')  ## error


geterrmessage()


## Import the log4r package.
library('log4r')

## Create a new logger object with create.logger().
logger <- create.logger()

## Set the logger's file output: currently only allows flat files.
logfile(logger) <- file.path('base.log')

## Set the current level of the logger.
level(logger) <- "INFO"

## Try logging messages at different priority levels.
debug(logger, 'A Debugging Message') ## Won't print anything
info(logger, 'An Info Message')
warn(logger, 'A Warning Message')
error(logger, 'An Error Message')
fatal(logger, 'A Fatal Error Message')


n <- 10
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
    setTxtProgressBar(pb, i)
    Sys.sleep(0.5)
}
close(pb)


library("progress")
pb <- progress_bar$new(total = n)
for (i in 1:n) {
    pb$tick()
    Sys.sleep(0.5)
}


if (!condition) stop(...)


stopifnot(TRUE)
stopifnot(TRUE, FALSE)


f <- function(x) {
    stopifnot(is.numeric(x), length(x) == 1)
    invisible(TRUE)
}

f(1)
f("1")
f(1:2)
f(letters)


x <- "1"
library("assertthat")
stopifnot(is.numeric(x))
assert_that(is.numeric(x))
assert_that(length(x) == 2)


a <- sqrt(2)
a * a == 2
a * a - 2


1L + 2L == 3L
1.0 + 2.0 == 3.0
0.1 + 0.2 == 0.3


all.equal(0.1 + 0.2, 0.3)
all.equal(0.1 + 0.2, 3.0)
isTRUE(all.equal(0.1 + 0.2, 3)) ## when you just want TRUE/FALSE


1 == NULL
all.equal(1, NULL)
identical(1, NULL)
identical(1, 1.)   ## TRUE in R (both are stored as doubles)
all.equal(1, 1L)
identical(1, 1L)   ## stored as different types


col_means <- function(df) {
  numeric <- sapply(df, is.numeric)
  numeric_cols <- df[, numeric]
  data.frame(lapply(numeric_cols, mean))
}


col_means(mtcars)
col_means(mtcars[, 0])
col_means(mtcars[0, ])
col_means(mtcars[, "mpg", drop = FALSE])
col_means(1:10)
col_means(as.matrix(mtcars))
col_means(as.list(mtcars))

mtcars2 <- mtcars
mtcars2[-1] <- lapply(mtcars2[-1], as.character)
col_means(mtcars2)


e <- function(i) {
  x <- 1:4
  if (i < 5) x[1:2]
  else x[-1:2]
}
f <- function() sapply(1:10, e)
g <- function() f()


g()
traceback()


debug(g)
g()


e <- function(i) {
  x <- 1:4
  if (i < 5) x[1:2]
  else x[-1:2] # oops! x[-(1:2)]
}
f <- function() sapply(1:10, e)
g <- function() f()


## make sure you have the 'sequences' package.
## Get readFasta2, the function to debug
library(devtools)
install_github("lgatto/sequences")
sequences:::debugme()
## Get an example file
f <- dir(system.file("extdata", package = "sequences"),
         full.names=TRUE, pattern = "moreDnaSeqs.fasta")
## BANG!
readFasta2(f)


f <- function() {
    x <- "1"
    log(x)
    message("x was the ", class(x), " ", x)
}
f()


f <- function() {
    x <- "1"
    try(log(x))
    message("x was the ", class(x), " ", x)
}
f()


try({
    a <- 1
    b <- "2"
    a + b
})


success <- try(1 + 2)
failure <- try(1 + "2", silent = TRUE)
class(success)
class(failure)


inherits(failure, "try-error")

if (inherits(failure, "try-error"))
	message("There was an error here.")


el <- list(1:10, c(-1, 1), TRUE, "1")
res <- lapply(el, log)
res
res <- lapply(el, function(x) try(log(x)))
res


default <- NULL
try(default <- read.csv("possibly-bad-input.csv"), silent = TRUE)


f <- function(x)
    if (x == 1) stop("Error!") else 1

f(1)
f(2)


safef <- failwith(NULL, f)
safef(1)
safef(2)


f <- function() {
    x <- "1"
    tryCatch(log(x),
             error = function(e) cat("There was an error!\n"))
    message("x was the ", class(x), " ", x)
}
f()


show_condition <- function(code) {
  tryCatch(code,
    error = function(c) "error",
    warning = function(c) "warning",
    message = function(c) "message"
  )
}
show_condition(stop("!"))
show_condition(warning("?!"))
show_condition(message("?"))
show_condition(0)


read.csv2 <- function(file, ...) {
  tryCatch(read.csv(file, ...), error = function(c) {
    c$message <- paste0(c$message, " (in ", file, ")")
    stop(c)
  })
}
read.csv("code/dummy.csv")
read.csv2("code/dummy.csv")


f <- function(x = 10) {
    lapply(seq_len(x), function(i) {
        ## make an example 2x2 contingency table
        d <- matrix(sample(4:10, 4), nrow = 2, ncol = 2)
        ## will produce warning if there is a 5 or less
        ## in the contingency table
        chisq.test(d)
    })
}


set.seed(1)
f()
set.seed(1)
withCallingHandlers(f(), warning=function(e) recover())


f <- function() g()
g <- function() h()
h <- function() stop("!")

tryCatch(f(), error = function(e) print(sys.calls()))
withCallingHandlers(f(), error = function(e) print(sys.calls()))


safelog <- function(x) {
  tryCatch(log(x),
           error = function(e) paste("an error with input", x),
           warning = function(e) paste("a warning with input", x))
}


log(1)
safelog(1)
log(-1)
safelog(-1)
log("a")
safelog("a")


## Report whenever e invoked
trace(sum)
hist(rnorm(100))
untrace(sum)


## Evaluate arbitrary code whenever e invoked
trace(e, quote(cat("i am", i, "\n")))
## Another way to enter browser whenver e invoked
trace(e, browser)
## stop tracing
untrace(e)


f <- function() {
    ## make an example 2x2 contingency table
    d <- matrix(sample(4:10, 4), nrow=2, ncol=2)
     chisq.test(d)
}
set.seed(1)
f() ## no warning

set.seed(11)
f() ## warning


if (any(d < 5))
  browser()


as.list(body(f))


trace("f", quote(if (any(d < 5)) browser()), at = 3)


f
body(f)


set.seed(1)
f() ## normal execution

set.seed(11)
f() ## enters browser mode


library("MSnbase")
data(itraqdata)
x <- itraqdata[[1]]
plot(x, full=TRUE)


debug(plot)
plot(x, full=TRUE)


trace("plot", browser,
      signature = c("Spectrum", "missing"))
plot(x, full=TRUE)

