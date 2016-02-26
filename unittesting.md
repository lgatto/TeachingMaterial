---
title: "Unit testing"
author: "Laurent Gatto"
---

These exercises were written by Martin Morgan and Laurent Gatto for a
[Bioconductor Developer Day workshop](http://bioconductor.org/help/course-materials/2013/BioC2013/developer-day-debug/).

# Introduction

> Whenever you are templted to type something into a print statement
> or a debugger expression, write it as a test insted -- Martin Fowler

**Why unit testing?**

- Writing code to test code;
- anticipate bugs, in particular for edge cases;
- anticipate disruptive updates;
- document and test observed bugs using specific tests.

Each section provides a function that supposedly works as expected,
but quickly proves to misbehave. The exercise aims at first writing
some dedicated testing functions that will identify the problems and
then update the function so that it passes the specific tests. This
practice is called unit testing and we use the RUnit package for
this.

See the
[Unit Testing How-To](http://bioconductor.org/developers/how-to/unitTesting-guidelines/)
guide for details on unit testing using the
[`RUnit`](http://cran.r-project.org/web/packages/RUnit/index.html)
package. The
[`testthat`](http://cran.r-project.org/web/packages/testthat/) is
another package that provides unit testing infrastructure. Both
packages can conveniently be used to automate unit testing within
package testing. 

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
## [1] "B" "Y" "A" "T" "J"
```
But


```r
## Bug!
isIn(c(x, "a"), LETTERS)
```

```
## [1] "B" "Y" "A" "T" "J" NA
```

### Solution

Write a unit test that demonstrates the issue


```r
## Unit test:
library("RUnit")
```

```
## Loading required package: methods
```

```r
test_isIn <- function() {
    x <- c("A", "B", "Z")
    checkIdentical(x, isIn(x, LETTERS))
    checkIdentical(x, isIn(c(x, "a"), LETTERS))

}

test_isIn()
```

```
## Error in checkIdentical(x, isIn(c(x, "a"), LETTERS)): FALSE 
## 
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

## The `testthat` syntax

`expect_that(object_or_expression, condition)` with conditions
- equals: `expect_that(1+2,equals(3))` or `expect_equal(1+2,3)`
- gives warning: `expect_that(warning("a")`, `gives_warning())`
- is a: `expect_that(1, is_a("numeric"))` or `expect_is(1,"numeric")`
- is true: `expect_that(2 == 2, is_true())` or `expect_true(2==2)`
- matches: `expect_that("Testing is fun", matches("fun"))` or `expect_match("Testing is fun", "f.n")`
- takes less: `than expect_that(Sys.sleep(1), takes_less_than(3))`

and

```r
test_that("description", {
    a <- foo()
    b <- bar()
    expect_equal(a, b)
})
```

## Batch unit testing

```r
library("testthat")
test_dir("./unittests/")
test_file("./unittests/test_foo.R")
```

# Exercises

## Column means

## Problem

The `col_means` function computes the means of all numeric columns in
a data frame (example from *Advanced R*, to illustrate defensive
programming).


```r
col_means <- function(df) {
  numeric <- sapply(df, is.numeric)
  numeric_cols <- df[, numeric]
  data.frame(lapply(numeric_cols, mean))
}

## Expected
col_means(mtcars)
```

```
##        mpg    cyl     disp       hp     drat      wt     qsec     vs
## 1 20.09062 6.1875 230.7219 146.6875 3.596563 3.21725 17.84875 0.4375
##        am   gear   carb
## 1 0.40625 3.6875 2.8125
```

```r
## Bugs
col_means(mtcars[, "mpg"])
```

```
## Error in df[, numeric]: incorrect number of dimensions
```

```r
col_means(mtcars[, "mpg", drop=FALSE])
```

```
##   X21 X21.1 X22.8 X21.4 X18.7 X18.1 X14.3 X24.4 X22.8.1 X19.2 X17.8 X16.4
## 1  21    21  22.8  21.4  18.7  18.1  14.3  24.4    22.8  19.2  17.8  16.4
##   X17.3 X15.2 X10.4 X10.4.1 X14.7 X32.4 X30.4 X33.9 X21.5 X15.5 X15.2.1
## 1  17.3  15.2  10.4    10.4  14.7  32.4  30.4  33.9  21.5  15.5    15.2
##   X13.3 X19.2.1 X27.3 X26 X30.4.1 X15.8 X19.7 X15 X21.4.1
## 1  13.3    19.2  27.3  26    30.4  15.8  19.7  15    21.4
```

```r
col_means(mtcars[, 0])
```

```
## Error in .subset(x, j): invalid subscript type 'list'
```

```r
col_means(mtcars[0, ])
```

```
##   mpg cyl disp  hp drat  wt qsec  vs  am gear carb
## 1 NaN NaN  NaN NaN  NaN NaN  NaN NaN NaN  NaN  NaN
```

```r
col_means(as.list(mtcars))
```

```
## Error in df[, numeric]: incorrect number of dimensions
```

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
## Warning in grep(x, y): argument 'pattern' has length > 1 and only the first
## element will be used
```

```
## [1] "abc" "a"
```

<!-- ### Solution -->



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

### Solution


```r
## Unit test:
library("RUnit")
test_ifcond <- function() {
    checkIdentical(5, ifcond(3, 2))
    checkIdentical(8, ifcond(2, 2))
    checkIdentical(5, ifcond(1, 2))
    checkIdentical(c(5, 8, 5), ifcond(3:1, c(2, 2, 2)))
}

test_ifcond()
```

```
## Warning in if (x > y) {: the condition has length > 1 and only the first
## element will be used
```

```
## Error in checkIdentical(c(5, 8, 5), ifcond(3:1, c(2, 2, 2))): FALSE 
## 
```

```r
## updated function:
ifcond <- function(x, y)
    ifelse(x > y, x*x - y*y, x*x + y*y)

test_ifcond()
```

```
## [1] TRUE
```

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
##              x           y
## [1,]  2.301199  0.60939751
## [2,]  1.380715  1.49984771
## [3,]  1.237685  0.58034293
## [4,] -1.723762 -0.45066731
## [5,] -1.389801 -0.05480788
```

```r
(p <- m[1, ])
```

```
##         x         y 
## 2.3011991 0.6093975
```

```r
distances(p, m)
```

```
## [1] 0.000000 1.280700 1.063910 4.162217 3.750287
```

```r
## Bug!
(dd <- data.frame(x, y))
```

```
##           x           y
## 1  2.301199  0.60939751
## 2  1.380715  1.49984771
## 3  1.237685  0.58034293
## 4 -1.723762 -0.45066731
## 5 -1.389801 -0.05480788
```

```r
(q <- dd[1, ])
```

```
##          x         y
## 1 2.301199 0.6093975
```

```r
distances(q, dd)
```

```
##   x
## 1 0
```

### Solution


```r
## Unit test:
library("RUnit")
test_distances <- function() {
    x <- y <- c(0, 1, 2)
    m <- cbind(x, y)
    p <- m[1, ]
    dd <- data.frame(x, y)
    q <- dd[1, ]
    expct <- c(0, sqrt(c(2, 8)))
    checkIdentical(expct, distances(p, m))
    checkIdentical(expct, distances(q, dd))
}

test_distances()
```

```
## Error in checkIdentical(expct, distances(q, dd)): FALSE 
## 
```

```r
## updated function
distances <- function(point, pointVec) {
    point <- as.numeric(point)
    x <- point[1]
    y <- point[2]
    xVec <- pointVec[,1]
    yVec <- pointVec[,2]
    dist <- sqrt((xVec - x)^2 + (yVec - y)^2)
    return(dist)
}

test_distances()
```

```
## [1] TRUE
```

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

### Solution


```r
## Unit test:
library(RUnit)
test_sqrtabs <- function() {
    checkIdentical(c(2, 0, 2), sqrtabs(c(-4, 0, 4)))
    checkIdentical(numeric(), sqrtabs(numeric()))
}
test_sqrtabs()
```

```
## Error in checkIdentical(numeric(), sqrtabs(numeric())): FALSE 
## 
```

```r
## updated function:
sqrtabs <- function(x) {
  v <- abs(x)
  sapply(seq_along(v), function(i) sqrt(v[i]))
}
test_sqrtabs()                          # nope!
```

```
## Error in checkIdentical(numeric(), sqrtabs(numeric())): FALSE 
## 
```

```r
sqrtabs <- function(x) {
  v <- abs(x)
  vapply(seq_along(v), function(i) sqrt(v[i]), 0)
}
test_sqrtabs()                          # yes!
```

```
## [1] TRUE
```

# Unit testing in a package 


## In a package

1. Create a directory `./mypackage/tests`.
2. Create the `testthat.R` file

```r
library("testthat")
library("mypackage")
test_check("sequences")
```

3. Create a sub-directory `./mypackage/tests/testthat` and include as
   many unit test files as desired that are named with the `test_`
   prefix and contain unit tests.

4. Suggest the unit testing package in your `DESCRIPTION` file:

```
Suggests: testthat
```

## Example from the `sequences` package

From the `./sequences/tests/testthat/test_sequences.R` file:

### Object creation and validity

We have a fasta file and the corresponding `DnaSeq` object.

1. Let's make sure that the `DnaSeq` instance is valid, as changes in
   the class definition might have altered its validity.

2. Let's verify that `readFasta` regenerates and identical `DnaSeq`
   object given the original fasta file.

```r
test_that("dnaseq validity", {
  data(dnaseq)
  expect_true(validObject(dnaseq))
})

test_that("readFasta", {
  ## loading _valid_ dnaseq
  data(dnaseq)
  ## reading fasta sequence
  f <- dir(system.file("extdata",package="sequences"),pattern="fasta",full.names=TRUE)
  xx <- readFasta(f[1])
  expect_true(all.equal(xx, dnaseq))
})
```

### Multiple implementations

Let's check that the R, C and C++ (via `Rcpp`) give the same result

```r
test_that("ccpp code", {
  gccountr <-
    function(x) tabulate(factor(strsplit(x, "")[[1]]))
  x <- "AACGACTACAGCATACTAC"
  expect_true(identical(gccount(x), gccountr(x)))
  expect_true(identical(gccount2(x), gccountr(x)))
})
```

## Exercise

Choose any data package of your choice and write a unit test that
tests the validity of all the its data. 

Tips 

- To get all the data distributed with a package, use `data(package = "packageName")`


```r
library("pRolocdata")
data(package = "pRolocdata")
```

- To test the validity of an object, use `validObject`




```r
data(andy2011)
validObject(andy2011)
```

```
## [1] TRUE
```

- Using the `testthat` syntax, the actual test for that data set would be 


```r
library("testthat")
expect_true(validObject(andy2011))
```


## Testing coverage in a package

The [covr](https://github.com/jimhester/covr) package:

![package coverage](./figs/covr.png)

We can use `type="all"` to examine the coverage in unit tests, examples and vignettes. This can
also be done interactively with Shiny:


```r
library(covr)
coverage <- package_coverage("/path/to/package/source", type="all")
shine(coverage)
```

[Coverage for all Bioconductor packages](https://codecov.io/github/Bioconductor-mirror).
