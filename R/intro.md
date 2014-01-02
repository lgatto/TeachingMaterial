# Introduction

## Pre-requisites

`variables` are used to store values: they bind a value to a name. The
content of a variable can be accessed by typing its name on the R
console prompt. Its content can also be overwritten.



```r
x
```

```
## Error: object 'x' not found
```

```r
x <- 1
y <- 1:10
z <- "foo"
x
```

```
## [1] 1
```

```r
x <- 2
x
```

```
## [1] 2
```


All these variables live in your works the global environment. They
can be listed with `ls()`. Variables can be deleted with
`rm()`. `ls()` and `rm()` are functions that take an input (as in
`rm(x)`; `ls()` is called without any explicit input values) and
either return a value (`ls()` returns the names of all existing
variables) or have side effects (`rm(x)` deletes the variable `x`).



```r
environment()
```

```
## <environment: R_GlobalEnv>
```

```r
ls()
```

```
## [1] "x" "y" "z"
```

```r
rm(x)
ls()
```

```
## [1] "y" "z"
```


`R` comes with a *limited* functionality. Thousands of third-pary
**packages** can be downloaded from on-line repositories (like CRAN)
and automatically installed with
`install.packages("newpackage")`. Packages are installed in specific
directories called **libraries** and can be loaded using
`library("newpackage")`.


```r
install.packages("devtools")
library("devtools")
```


To see what packages are available, one can use `installed.packages()`.


```r
allpkgs <- installed.packages()
nrow(allpkgs)
```

```
## [1] 841
```



```r
head(allpkgs)
```


`R` comes with a extensive help system. Type `helps.start()` to start
the HTML version of `R`'s online documentation. Each function comes
with its dedicated manual that can be accessed using `help("ls")` or
`?ls`. Reading the `ls` manual, we see that we can provide it with a
`list` input, i.e. *a character vector naming objects to be
removed*. Below, we generate this list of variable (object) names with
`ls()` and pass it directly as input to `rm()` to get a clean
workspace.


```r
ls()
```

```
## [1] "allpkgs" "y"       "z"
```

```r
rm(list = ls())
ls()
```

```
## character(0)
```


To get all the details about a package and its documentation:


```r
packageDescription("devtools")
help(package = "devtools)"
```


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
class(a)
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
class(b)
```

```
## [1] "numeric"
```

```r

a + b
```

```
## [1] 2
```

### a functional programming language 

Functions are *first class citizens*. Use function (outputs) as inputs
to other functions (see the `rm(list = ls())` example above), to create other variables, ... we will make use of this throughout the workshop. More details [here](https://github.com/lgatto/R-functional-programming#readme).

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

|        | dimensions | content types |
|--------|------------|---------------|
| vector | 1          | 1             |
| matrix | 2          | 1             |
| array  | n          | 1             |
| list   | 1          | `length(l)`   |
| data.frame | 2      | `ncol(dfr)`   |

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


