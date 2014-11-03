
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
## [1] "Z" "X" "T" "B" "G"
```
But


```r
## Bug!
isIn(c(x, "a"), LETTERS)
```

```
## [1] "Z" "X" "T" "B" "G" NA
```
