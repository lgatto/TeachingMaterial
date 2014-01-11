
## ----, echo=FALSE--------------------------------------------------------
opts_chunk$set(eval=FALSE)


## ------------------------------------------------------------------------
library("camweather")
f <- weatherfile("2014-01-01")
f
w <- read.table(f, header = FALSE,
                comment.char = "#",
                sep = "\t")
dim(w)
head(w)


## ------------------------------------------------------------------------
hd <- readLines(f)
hd <- hd[grep("#", hd)]
hd <- sub("#", "", hd)
hd <- hd[7:8]
hd <- gsub(" ", "", hd)
hd <- strsplit(hd, "\t")
hd <- paste0(hd[[1]], " [", hd[[2]], "]")
hd <- sub(" \\[\\]", "", hd)
names(w) <- hd


## ------------------------------------------------------------------------
class(w$Time)
w$Time <- strptime(paste(basename(f), w$Time), "%Y_%m_%d %H:%M")
w$Day <- as.Date(basename(f), "%Y_%m_%d")
class(w$Time)
summary(w)


## ----, echo=TRUE, eval=TRUE----------------------------------------------
library("camweather")
w <- weatherdata("2014-01-01")
w <- nounits(w)


## ----, eval=TRUE---------------------------------------------------------
par(mfrow = c(2, 2))
plot(w$Time, w[, "Temp"], type = "b", xlab = "Time", ylab = "Temp")
plot(w$Time, w[, "WindSp"], type = "b", xlab = "Time", ylab = "Wind speed")
plot(w$Time, w[, "Rain"], type = "b", xlab = "Time", ylab = "Rain")
plot(w$Time, w[, "Press"], type = "b", xlab = "Time", ylab = "Pressure")


## ----, eval = TRUE, tidy = FALSE-----------------------------------------
boxplot(w$WindSp ~ w$WindDr) ## NOT boxplot(w$WindDr, w[, "WindSp"])
pairs(w[, c(2, 5, 6, 9)])


## ----, eval = TRUE-------------------------------------------------------
temp0 <- w[, "Temp"]
temp <- temp0 - min(temp0) ## min is 0
temp <- temp/max(temp) ## max is 1

press0 <- w[, "Press"]
press <- press0 - min(press0)
press <- press/max(press)


## ----, eval=TRUE---------------------------------------------------------
## Plot with minimal decoration
par(mar = c(5, 4, 2, 4))
plot(w$Time,  temp , type = "l",
     xlab = "Time", ylab = "Temp [deg C]",
     yaxt = "n", col = "steelblue")
lines(w$Time, press, col = "red")

## Axis, title and legends
axis(2, at = seq(0, 1, length = 11),
     labels = seq(min(temp0),
         max(temp0),
         length = 11))
axis(4, at = seq(0, 1, length = 11),
     labels = seq(min(press0),
         max(press0),
         length = 11))
mtext("Pressure [mBar]", 4, line = 3)
title("2014-01-01")
legend("top", c("Temperature", "Pressure"),
       col = c("steelblue", "red"), lty = 1,
       bty = "n")


## ----, eval=FALSE--------------------------------------------------------
## write.csv(w, file = "w.csv")


## ----, eval = TRUE-------------------------------------------------------
x <- sqrt(2)
x^2 == 2


## ----, eval = TRUE-------------------------------------------------------
all.equal(x^2, 2)


## ------------------------------------------------------------------------
choices <- c("RStudio", "Wordpad", "emacs", "vim", "Notepad++")
mychoice <- menu(choices, graphics = FALSE, title = "Best editor ever")
cat("Best editor ever is ", choices[mychoice], "\n")


## ------------------------------------------------------------------------
for (i in 1:3)
    print(i + 1)

k <- 3
while (k > 0) {
    print(k^2)
    k <- k - 1
}


## ----, eval = TRUE-------------------------------------------------------
x <- numeric()
n <- length(x)

for (i in 1:n)
    print(i)

for (i in seq_len(n))
    print(i)


## ------------------------------------------------------------------------
k <- 3
repeat {
    print(k)
    if (k == 0) break
    k <- k - 1
}


## ----, eval=TRUE---------------------------------------------------------
x <- 1:5
ifelse(x < 3, 1, 2)
ifelse(x < 3, 10:15, 20:25)
ifelse(x < 3, 10:15, -1)


## ----, eval=TRUE---------------------------------------------------------
m <- matrix(1:6, ncol = 2)
m
ifelse(m < 2, m, c(10, 11, 12))


## ------------------------------------------------------------------------
switch(1+0, 1, 2, 3)
switch(1+1, 1, 2, 3)


## ------------------------------------------------------------------------
switch(letters[1+0], a = 1, b = 2, c = 3)
switch(letters[1+2], a = 1, b = 2, c = 3)


## ------------------------------------------------------------------------
l <- list(1:4, letters, month.name)
sapply(l, length)


## ------------------------------------------------------------------------
lapply(l, length)


## ------------------------------------------------------------------------
lapply(c(3, 4, 5), seq_len)
sapply(c(3, 4, 5), seq_len)


## ------------------------------------------------------------------------
m <- matrix(rnorm(30), nrow = 6)
apply(m, 1, median)
apply(m, 2, median)


## ------------------------------------------------------------------------
mapply(rep, 1:4, 4:1)


## ------------------------------------------------------------------------
x <- 1:12
k <- rep(letters[1:3], 4)
tapply(x, k, sum)


## ------------------------------------------------------------------------
myfun <- function(x, y) {
    a <- x^2
    b <- sqrt(y)
    res <- a/b
    return(res)
}

myfun(2, 1)
myfun(4, 2)


## ------------------------------------------------------------------------
myfun2 <- function(x, y = 2) {
    a <- x^2
    b <- sqrt(y)
    res <- a/b
    return(res)
}

myfun2(4, 2)
myfun2(4)


## ------------------------------------------------------------------------
myfun(c(4, 4), c(1, 2))


## ------------------------------------------------------------------------
centre <- function(x, type) {
       switch(type,
              mean = mean(x),
              median = median(x),
              trimmed = mean(x, trim = .1))
     }
x <- rcauchy(10)
centre(x, "mean")
centre(x, "median")
centre(x, "trimmed")


## ------------------------------------------------------------------------
plot1toN <- function(n, ...) plot(1, n, ...)


## ------------------------------------------------------------------------
m <- matrix(rnorm(12), ncol = 4)
apply(m, 1, function(x) sum(x^2))


## ----, eval=TRUE---------------------------------------------------------
x <- 1
f <- function(x) {
    x <- x + 1
    return(x)
}
f(x)
x ## unchanged


## ----, eval = TRUE-------------------------------------------------------
g <- function() {
    x <- x + 1
    return(x)
}
x <- 1
g()
x ## unchanged


## ----, eval=TRUE---------------------------------------------------------
rm(x)
g()


## ------------------------------------------------------------------------
myenv <- new.env()
myenv$x <- 1
updatex <- function(e, newx)
    assign("x", newx, envir = e)
myenv$x
updatex(myenv, 10)
myenv$x


## ------------------------------------------------------------------------
##' Get the weather data for a day.
##'
##' Data are immediate at \code{Time} except wind speed (average since
##' previous \code{Time}) and wind direction (most frequent since
##' previous \code{Time}.)  Sun and rain values are cumulative from
##' code{Start}. \code{MxWSpd} gives max wind speed since previous
##' \code{Time}.
##' 
##' @title Weather data
##' @param date A character describing a date with format
##' \code{"YYYY-MM-DD"}.
##' @return A \code{data.frame} with the weather data for the
##' corresponding \code{date}.
##' @author Laurent Gatto <lg390@@cam.ac.uk>
##' @seealso \code{\link{nounits}} to remove the units from the
##' \code{data.frame}'s names.
##' @examples
##' x <- weatherdata("2012-12-25")
##' dim(x)
##' head(x)
##' plot(x$Time, x[, "Temp [degC]"], type = "b")
weatherdata <- function(date) {
    f <- weatherfile(date)
    if (length(f) > 1) {
        warning("Found ", length(f), " files. Using first one ",
                basename(f[1]))
        f <- f[1]
    }
    w <- read.table(f, header = FALSE,
                    comment.char = "#",
                    sep = "\t")
    hd <- readLines(f)
    hd <- hd[grep("#", hd)]
    hd <- sub("#", "", hd)
    hd <- hd[7:8]
    hd <- gsub(" ", "", hd)
    hd <- strsplit(hd, "\t")
    hd <- paste0(hd[[1]], " [", hd[[2]], "]")
    hd <- sub(" \\[\\]", "", hd)
    names(w) <- hd
    w$Time <- strptime(paste(basename(f), w$Time), "%Y_%m_%d %H:%M")
    w$Day <- as.Date(basename(f), "%Y_%m_%d")
    return(w)
}


## ----, eval = TRUE-------------------------------------------------------
X <- rnorm(1e6)
f <- function(x, k = .8) mean(x, trim = k)
f(X)
system.time(f(X))
summary(replicate(10, system.time(f(X))["elapsed"]))


## ----, eval=TRUE---------------------------------------------------------
n <- 1e4
f1 <- function(n) {
    l <- list()
    for (i in seq_len(n))
        l[[i]] <- seq(i)
    return(l)
}

f2 <- function(n) {
    l <- vector("list", length = n)
    for (i in seq_len(n))
        l[[i]] <- seq(i)
    return(l)
}

f3 <- function(n) 
    lapply(seq_len(n), seq)


## ----, tidy=FALSE, eval=TRUE, cache=TRUE---------------------------------
library("rbenchmark")
benchmark(f1(n), f2(n), f3(n),
          columns = c("test", "replications", "elapsed", "relative"),
          replications = 10)


