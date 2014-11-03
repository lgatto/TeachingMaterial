
# Introduction

Why unit testing?

Each section provides a function that supposedly works as expected,
but quickly proves to misbehave. The exercise aims at first writing
some dedicated testing functions that will identify the problems and
then update the function so that it passes the specific tests. This
practice is called unit testing and we use the RUnit package for
this. See the
[Unit Testing How-To](http://bioconductor.org/developers/how-to/unitTesting-guidelines/)
guide for details on unit testing using RUnit.

# Exercises

## Subsetting

### Problem

This function should return the elements of `x` that are in `y`.


```r
## Example
isIn <- function(x, y) {
    sel <- match(x, y)
    y[sel]
}

## Expected
x <- sample(LETTERS, 5)
isIn(x, LETTERS)
```

```
## [1] "X" "M" "H" "Y" "O"
```
But


```r
## Bug!
isIn(c(x, "a"), LETTERS)
```

```
## [1] "X" "M" "H" "Y" "O" NA
```

### Solution


```r
## Unit test:
library("RUnit")
test_isIn <- function() {
    x <- c("A", "B", "Z")
    checkIdentical(x, isIn(x, LETTERS))
    checkIdentical(x, isIn(c(x, "a"), LETTERS))

}

test_isIn()
```

```
## Error in checkIdentical(x, isIn(c(x, "a"), LETTERS)): FALSE
```

Unpdate the buggy function until the unit test succeeds


```r
## updated function
isIn <- function(x, y) {
    sel <- x %in% y
    x[sel]
}

test_isIn()
```

```
## [1] TRUE
```


## Character matching

### Problem


```r
isExactIn <- function(x, y)
    y[grep(x, y)]

## Expected
isExactIn("a", letters)
```

```
## [1] "a"
```

```r
## Bugs
isExactIn("a", c("abc", letters))
```

```
## [1] "abc" "a"
```

```r
isExactIn(c("a", "z"), c("abc", letters))
```

```
## Warning in grep(x, y): argument 'pattern' has length > 1 and only the
## first element will be used
```

```
## [1] "abc" "a"
```

## If conditions with length > 1

### Problem


```r
ifcond <- function(x, y) {
    if (x > y) {
        ans <- x*x - y*y
    } else {
        ans <- x*x + y*y
    } 
    ans
}

## Expected
do(3, 2)
```

```
## Error in eval(expr, envir, enclos): could not find function "do"
```

```r
do(2, 2)
```

```
## Error in eval(expr, envir, enclos): could not find function "do"
```

```r
do(1, 2)
```

```
## Error in eval(expr, envir, enclos): could not find function "do"
```

```r
## Bug!
do(3:1, c(2, 2, 2))
```

```
## Error in eval(expr, envir, enclos): could not find function "do"
```

## Know your inputs

### Problem


```r
## Example
distances <- function(point, pointVec) {
    x <- point[1]
    y <- point[2]
    xVec <- pointVec[,1]
    yVec <- pointVec[,2]
    sqrt((xVec - x)^2 + (yVec - y)^2)
}

## Expected
x <- rnorm(5)
y <- rnorm(5)

m <- cbind(x, y)
p <- m[1, ]

distances(p, m)
```

```
## [1] 0.0000000 3.5610160 0.9723341 2.0717517 2.6811202
```

```r
## Bug!
dd <- data.frame(x, y)
q <- dd[1, ]

distances(q, dd)
```

```
##   x
## 1 0
```

## Iterate on 0 length

### Problem


```r
sqrtabs <- function(x) {
    v <- abs(x)
    sapply(1:length(v), function(i) sqrt(v[i]))
}

## Expected
all(sqrtabs(c(-4, 0, 4)) == c(2, 0, 2))
```

```
## [1] TRUE
```

```r
## Bug!
sqrtabs(numeric())
```

```
## [[1]]
## [1] NA
## 
## [[2]]
## numeric(0)
```
