# Introduction

## Pre-requisites


## `R` is

### a dynamic language


```r
a <- 1L
mode(a)
```

```
## [1] "numeric"
```

```r
typeof(a)
```

```
## [1] "integer"
```

```r

b <- 1
mode(b)
```

```
## [1] "numeric"
```

```r
typeof(b)
```

```
## [1] "double"
```

```r

x <- "1"
mode(x)
```

```
## [1] "character"
```

```r
typeof(x)
```

```
## [1] "character"
```

```r

paste(a, x)
```

```
## [1] "1 1"
```

### a functional programming language 

Functions are *first class citizens*. Use function (outputs) as inputs
to other variables, to create other variables, ... we will make use of
this throughout the workshop. More details
[here](https://github.com/lgatto/R-functional-programming#readme).

### lazy

Calling a function with argument *wait for 3 seconds* takes no time to return

```r
f <- function(x) {
    10
}
system.time(f(Sys.sleep(3)))
```

```
##    user  system elapsed 
##       0       0       0
```


unless we *force* evaluation


```r
f <- function(x) {
    force(x)
    10
}
system.time(f(Sys.sleep(3)))
```

```
##    user  system elapsed 
##   0.000   0.000   3.003
```


(Example originally from Hadley Wickham's devtools page)

### object-oriented

with multiple frameworks: S3, S4 and S4 reference classes

## Data structures

Introduce data structures and manipulation 

- vectors
- numerics, characters, logicals, factors
- matrix
- array
- list
- data.frame
- special values: `NULL`, `NA`, `NaN`

### Objects

As in OO programming:


```r
class(x <- rnorm(100))
```

```
## [1] "numeric"
```

```r
y <- rnorm(100)
class(1L)
```

```
## [1] "integer"
```

```r
class("123")
```

```
## [1] "character"
```

```r
class(sum)
```

```
## [1] "function"
```

```r
class(lm(y ~ x))
```

```
## [1] "lm"
```


