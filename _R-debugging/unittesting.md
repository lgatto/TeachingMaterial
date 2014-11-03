These exercises were written by Martin Morgan and Laurent Gatto to be
used during
[Bioconductor Developer Day workshop](http://bioconductor.org/help/course-materials/2013/BioC2013/developer-day-debug/)
activities.

# Introduction

**Why unit testing?**

Each section provides a function that supposedly works as expected,
but quickly proves to misbehave. The exercise aims at first writing
some dedicated testing functions that will identify the problems and
then update the function so that it passes the specific tests. This
practice is called unit testing and we use the RUnit package for
this. See the
[Unit Testing How-To](http://bioconductor.org/developers/how-to/unitTesting-guidelines/)
guide for details on unit testing using RUnit.

# Example

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
## [1] "S" "T" "P" "M" "C"
```
But


```r
## Bug!
isIn(c(x, "a"), LETTERS)
```

```
## [1] "S" "T" "P" "M" "C" NA
```

### Solution

Write a unit test that demonstrates the issue


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

Update the buggy function until the unit test succeeds


```r
## updated function
isIn <- function(x, y) {
    sel <- x %in% y
    x[sel]
}

test_isIn() ## the bug is fixed and monitored
```

```
## [1] TRUE
```

# Exercises

## Character matching

### Problem

What are the exact matches of `x` in `y`?


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

<!-- ### Solution -->

<!-- ```{r} -->
<!-- ## Unit test: -->
<!-- library("RUnit") -->
<!-- test_isExactIn <- function() { -->
<!--     checkIdentical("a", isExactIn("a", letters)) -->
<!--     checkIdentical("a", isExactIn("a", c("abc", letters))) -->
<!--     checkIdentical(c("a", "z"), isExactIn(c("a", "z"), c("abc", letters))) -->
<!-- } -->

<!-- test_isExactIn() -->

<!-- ## updated function: -->
<!-- isExactIn <- function(x, y) -->
<!--     x[x %in% y] -->

<!-- test_isExactIn() -->
<!-- ``` -->

## If conditions with length > 1

### Problem

If `x` is greater than `y`, we want the difference of their
squares. Otherwise, we want the sum.


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
ifcond(3, 2)
```

```
## [1] 5
```

```r
ifcond(2, 2)
```

```
## [1] 8
```

```r
ifcond(1, 2)
```

```
## [1] 5
```

```r
## Bug!
ifcond(3:1, c(2, 2, 2))
```

```
## Warning in if (x > y) {: the condition has length > 1 and only the first
## element will be used
```

```
## [1]  5  0 -3
```

<!-- ### Solution -->

<!-- ```{r} -->
<!-- ## Unit test: -->
<!-- library("RUnit") -->
<!-- test_ifcond <- function() { -->
<!--     checkIdentical(5, ifcond(3, 2)) -->
<!--     checkIdentical(8, ifcond(2, 2)) -->
<!--     checkIdentical(5, ifcond(1, 2)) -->
<!--     checkIdentical(c(5, 8, 5), ifcond(3:1, c(2, 2, 2))) -->
<!-- } -->

<!-- test_ifcond() -->

<!-- ## updated function: -->
<!-- ifcond <- function(x, y) -->
<!--     ifelse(x > y, x*x - y*y, x*x + y*y) -->

<!-- test_ifcond() -->
<!-- ``` -->

## Know your inputs

### Problem

Calculate the euclidean distance between a single point and a set of
other points.


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

(m <- cbind(x, y))
```

```
##               x          y
## [1,]  0.9337659 -1.7568430
## [2,] -0.2461883 -0.3463200
## [3,] -0.3517330 -0.4897523
## [4,] -1.0087053 -0.5456311
## [5,]  0.3394810  1.2328814
```

```r
(p <- m[1, ])
```

```
##          x          y 
##  0.9337659 -1.7568430
```

```r
distances(p, m)
```

```
## [1] 0.000000 1.838985 1.805000 2.289155 3.048217
```

```r
## Bug!
(dd <- data.frame(x, y))
```

```
##            x          y
## 1  0.9337659 -1.7568430
## 2 -0.2461883 -0.3463200
## 3 -0.3517330 -0.4897523
## 4 -1.0087053 -0.5456311
## 5  0.3394810  1.2328814
```

```r
(q <- dd[1, ])
```

```
##           x         y
## 1 0.9337659 -1.756843
```

```r
distances(q, dd)
```

```
##   x
## 1 0
```

<!-- ### Solution -->

<!-- ```{r} -->
<!-- ## Unit test: -->
<!-- library("RUnit") -->
<!-- test_distances <- function() { -->
<!--     x <- y <- c(0, 1, 2) -->
<!--     m <- cbind(x, y) -->
<!--     p <- m[1, ] -->
<!--     dd <- data.frame(x, y) -->
<!--     q <- dd[1, ] -->
<!--     expct <- c(0, sqrt(c(2, 8))) -->
<!--     checkIdentical(expct, distances(p, m)) -->
<!--     checkIdentical(expct, distances(q, dd)) -->
<!-- } -->

<!-- test_distances() -->

<!-- ## updated function -->
<!-- distances <- function(point, pointVec) { -->
<!--     point <- as.numeric(point) -->
<!--     x <- point[1] -->
<!--     y <- point[2] -->
<!--     xVec <- pointVec[,1] -->
<!--     yVec <- pointVec[,2] -->
<!--     dist <- sqrt((xVec - x)^2 + (yVec - y)^2) -->
<!--     return(dist) -->
<!-- } -->

<!-- test_distances() -->
<!-- ``` -->

## Iterate on 0 length

### Problem

Calculate the square root of the absolute value of a set of numbers.


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

<!-- ### Solution -->

<!-- ```{r} -->
<!-- ## Unit test: -->
<!-- library(RUnit) -->
<!-- test_sqrtabs <- function() { -->
<!--     checkIdentical(c(2, 0, 2), sqrtabs(c(-4, 0, 4))) -->
<!--     checkIdentical(numeric(), sqrtabs(numeric())) -->
<!-- } -->
<!-- test_sqrtabs() -->

<!-- ## updated function: -->
<!-- sqrtabs <- function(x) { -->
<!--   v <- abs(x) -->
<!--   sapply(seq_along(v), function(i) sqrt(v[i])) -->
<!-- } -->
<!-- test_sqrtabs()                          # nope! -->

<!-- sqrtabs <- function(x) { -->
<!--   v <- abs(x) -->
<!--   vapply(seq_along(v), function(i) sqrt(v[i]), 0) -->
<!-- } -->
<!-- test_sqrtabs()                          # yes! -->
<!-- ``` -->

