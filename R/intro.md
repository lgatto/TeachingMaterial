Introduction
===




# Pre-requisites

**Variables** are used to store values: they bind a value to a name. The
content of a variable can be accessed by typing its name on the R
console prompt. Its content can also be overwritten.



```r
x
x <- 1
y <- 1:10
z <- "foo"
x
x <- 2
x
```


All these variables live in your works the global environment. They
can be listed with `ls()`. Variables can be deleted with
`rm()`. `ls()` and `rm()` are functions that take an input (as in
`rm(x)`; `ls()` is called without any explicit input values) and
either return a value (`ls()` returns the names of all existing
variables) or have side effects (`rm(x)` deletes the variable `x`).



```r
environment()
ls()
rm(x)
ls()
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
rm(list = ls())
ls()
```


To get all the details about a package and its documentation:


```r
packageDescription("devtools")
help(package = "devtools)"
```


## Comments

Using the `#` character.

## Session information


```r
sessionInfo()
version
```


# `R` is

## a dynamic language


```r
a <- 1L
mode(a)
typeof(a)
class(a)

b <- 1
mode(b)
typeof(b)
class(b)

a + b
```

## a functional programming language 

Functions are *first class citizens*. Use function (outputs) as inputs
to other functions (see the `rm(list = ls())` example above), to
create other variables, ... we will make use of this throughout the
workshop. More details
[here](https://github.com/lgatto/R-functional-programming#readme).

## lazy

Calling a function with argument *wait for 3 seconds* takes no time to return


```r
f <- function(x) {
    10
}
system.time(f(Sys.sleep(3)))
```


unless we *force* evaluation


```r
f <- function(x) {
    force(x)
    10
}
system.time(f(Sys.sleep(3)))
```


(Example originally from Hadley Wickham's devtools page)

## object-oriented

with multiple frameworks: S3, S4 and S4 reference classes

# Data structures

|        | dimensions/length | content types |
|--------|------------|---------------|
| vector | 1          | 1             |
| matrix | 2          | 1             |
| array  | n          | 1             |
| list   | 1          | `length(l)`   |
| data.frame | 2      | `ncol(dfr)`   |

## Vectors

The basic type in `R` is a vector. Vectors have a specific length
(queried with the `length()` function) and have optional names
(accessed with `names()` or set with `names() <- `). Subsets of
vectors can be extracted by numerical indices (starting at 1) or names
(if available) using the `[` accessor.


```r
x <- 1
length(x)
x <- c(1, 2, 4)  ## combine values into a vector
length(x)
x[2]  ## second element of x
```



```r
names(x)
names(x) <- c("a", "b", "c")
x["b"]
```


We can specifically update subsets of a vector...


```r
x[2] <- x[2] * 10
x
```


... or delete specific elements using negative indices.


```r
y <- x[-2]
y
```


## Numerics, characters, logicals, factors

Vectors are of one unique type:


```r
typeof(c(1, 10, 1))
typeof(c(1L, 10L, 1L))
typeof(c(1L, 10, 1))
typeof(letters)
typeof(c(TRUE, FALSE))  ## TRUE and FALSE are reserved words
gender <- factor(c("male", "female", "male"))
gender
class(gender)
typeof(gender)
```

## Vectorised operations


```r
x <- c(1, 2, 3)
y <- c(3, 2, 1)
x
y
x + y
x^2
sqrt(x)
```



```r
x <- c("a", "b", "c")
paste(x, 1, sep = ".")
```


### Recycling


```r
x <- c(1, 2, 3, 4)
y <- c(1, 2)
x + y
z <- c(1, 2, 3)
x + z
```

## Generating vectors

- `seq` and **argument matching**
- `:`

### Exercise:

Using `rep`, how to generate 
- 3 repetitions of 1 to 10: 1, 2, ..., 9, 10, 1, 2, ..., 9, 10, 1, 2, ..., 9, 10
- repeating numbers 1 to 10 each 3 times: 1, 1, 1, 2, 2, 2, ..., 9, 9, 9, 10, 10, 10

## Matrix

A vector with 2 dimensions


```r
m <- c(1, 2, 3, 4, 5, 6)
m
dim(m) <- c(2, 3)  ## always rows first
m
class(m)  ## a matrix
mode(m)  ## of numerics
```


Accessing/subsetting is the same as vectors, keeping in mind that we have 2 dimensions:


```r
m[1:2, 1:2]
m[, 2:3]
m[1, 1] <- 10
m
m[, -1]
```


Naming matrices' rows and columns


```r
colnames(m) <- letters[1:3]
rownames(m) <- LETTERS[1:2]
m
m["A", "c"]
```


## Simplification/dropping of dimensions


```r
m[-1, ]  ## becomes a vector
m[-1, , drop = FALSE]  ## remains a matrix
m[, -(1:2)]  ## becomes a vector
m[, -(1:2), drop = FALSE]  ## remains a matrix
```


## Array

Like a matrix with > 2 dimensions.

## List

A `list` is a generic vector that can store elements of different
types. Lists can also have names (same syntax as vectors) and accessed
with `[[]]` or `$element.name`


```r
l <- list(M = matrix(1:10, ncol = 5), V = 1:2, F = sum)
l
l$M + l[["V"]]
l$F(l$M)  ## same as sum(1:10)
```


### Difference between `[` and ``[``

The former returns a sub-set of the list, the latter an individual elements


```r
class(l[2])
class(l[[2]])
```


Hence


```r
class(l[2:3])
class(l[[2:3]])
```


## Data.frame

A `data.frame` is a list whose elements are all of the same lengths
and that is represented as a table. Matrix-like subsetting `[,]` using
`names` (by definition `colnames`) and `rownames` or indices can be
used.

### Exercise

Generate a `data.frame` of patient data, including their first names,
surnames, age, gender, weights and whether they give consent for their
data to be made public.


```r
age <- c(50, 21, 35, 45, 28,
         31, 42, 33, 57, 62)
weight <- c(70.8, 67.9, 75.3, 61.9, 72.4,
            69.9, 63.5, 71.5, 73.2, 64.8)
firstName <- c("Adam", "Eve", "John", "Mary",
               "Peter", "Paul", "Joanna", "Matthew",
               "David", "Sally")
secondName <- c("Jones", "Parker", "Evans",
                "Davis", "Baker", "Daniels",
                "Edwards", "Smith", "Roberts", "Wilson")
consent <- c(TRUE, TRUE, FALSE, TRUE, FALSE,
             FALSE, FALSE, TRUE, FALSE, TRUE)
gender <- c("Male", "Female", "Male", "Female",
            "Male", "Male", "Female", "Male",
            "Male", "Female")
patients <- data.frame(firstName, secondName,
                       gender = factor(gender),
                       age, weight, consent,
                       stringsAsFactors=FALSE)
```


By default, `character`s are converted to `factor`s when generating
`data.frame`s; use `stringsAsFactors = FALSE` to keep the `characters`
as is and convert to `factor`s explicitly.

## Special values: `NULL`, `NA`, `NaN`


```r
class(NA)
class(NaN)
class(NULL)
```



```r
l[[2]] <- NULL
l
```



```r
sum(1, 2, NA)
sum(1, 2, NA, na.rm = TRUE)
```


## Coercing

The convert from one type to another, use `as.*`. To test whether a
variable is of a certain type, use `is*`. Type `as.TAB` and `is.TAB`
to tab-complete as list of functions.


```r
as.numeric("1")
as.character(1 + 2)
as.integer(1.9)  ## see also floor and ceiling
as.numeric("1a")
```


## Objects

As in OO programming:


```r
class(x <- rnorm(100))
y <- rnorm(100)
class(1L)
class("123")
class(sum)
class(lm(y ~ x))
```



```r
library("affydata")
data("Dilution")
class(Dilution)
```

 
