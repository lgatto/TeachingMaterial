Introduction
===




# Pre-requisites

## Basic tools

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


**Functions** are used to perform actions. They can be recognised by
the parenthesis following the function's name. Functions have optional
inputs, termed arguments, that are passed to the function through by
specifying then inside the parentheses. Functions also return output
values (that are generally computed based on the input values) and
sometimes also side effects, when they produce something else that
returning a value.


```r
sum(1, 2, 3)
sqrt(2)
x <- 2
x^2
```


## Workspace

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


## Working directory

`R` is running in a directory, called the working directory. That's
where `R` expects to find and save files. The working directory can be
queried with the `getwd()` function. It does not take any inputs
(there is nothing inside the `()` and returns the name the current
working directory).


```r
getwd()
```


The working directory can also be updated with
`setwd("path/to/workingdir")`, where `"path/to/workingdir"` points to
the new working directory. The working directory is typically your
current project dorectory, where input data, source code, result
files, figures, ... are stored.

It is of course possible to access and store files in other
directories by specifying the full path. 

## Packages

`R` comes with a *limited* functionality. Thousands of third-party
**packages** can be downloaded from dedicated on-line repositories
(like CRAN) and automatically installed with
`install.packages("newpackage")`. Packages are installed in specific
directories called **libraries** and can be loaded using
`library("newpackage")`.


```r
install.packages("devtools")
library("devtools")
```


It is of course possible to download and install packages using
`install.packages(, repo = NULL)`. However, this will force you to
handle all dependencies (and there can be many) manually too.


It is also possible to install packages from
[github](https://github.com/) using functionality from the package we
have just installed. We first need to load the package with the
`library` function to be able to use the functions that it
provides. Below, we install the `camweather` package (that provides
data about the weather in Cambridge) from the `lgatto` github
repository.


```r
library("devtools")
install_github("lgatto/camweather")
```


To see what packages are available, one can use `installed.packages()`.


```r
allpkgs <- installed.packages()
nrow(allpkgs)
```



```r
head(allpkgs)
```


## Help

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


Some packages provide an additional piece of documentation called a
`vignette`. It is generally displayed as a `pdf` file or, more
recently, as `html` pages. Vignettes are written as `LaTeX` (or
`markdown`) documents with interleaved code chunks. These are
evaluated and replaced by their output. These documents generally
present a high-level overview of the package's functionality, as
opposed to individual manual pages, that provide detailed information
about a single object.

To find out if a package has vignettes, use `vignette(package =
"packagename")`. To load a specific vignette, use
`vignette("vignettename", package = "packagename")` or
`vignette("vignettename")`.

## Comments

Using the `#` character.

## Session information


```r
sessionInfo()
version
```


# `R` is

## A dynamic language


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

## A functional programming language 

Functions are *first class citizens*. Use function (outputs) as inputs
to other functions (see the `rm(list = ls())` example above), to
create other variables, ... we will make use of this throughout the
workshop. More details
[here](https://github.com/lgatto/R-functional-programming#readme).

## Lazy

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

## Object-oriented

With multiple frameworks: S3, S4 and S4 reference classes.

# Data structures

|          | dimensions  | content types |
|----------|-------------|---------------|
| `vector` | 1: `lenght` | 1             |
| `matrix` | 2: `dim`    | 1             |
| `array`  | n: `dim`    | 1             |
| `list`   | 1: `length` | `length()`    |
| `data.frame` | 2: `dim`| `ncol()`      |

Is the data a collection of scalars of same *type* (defined later)?
 - yes:
   - linear?
     - yes `vector`
	 - no `matrix` or `array` (if > 2 dimensions)
 - no: `list`
   - elements of same length (and generally `vectors`): `data.frame`

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


## Vector modes: numerics, character, logical, factor

Vectors are of one unique type:

          | mode | class  | typeof |
----------|------|--------|--------| 
character | character | character | character |
logical   | logical | logical | logical |
numeric   | numeric | integer or numeric | integer of double | 
factor    | numeric | factor | integer | 



```r
typeof(c(1, 10, 1))
typeof(c(1L, 10L, 1L))
typeof(c(1L, 10, 1))
typeof(letters)
typeof(c(TRUE, FALSE))  ## TRUE and FALSE are reserved words
```

## Factors


```r
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
- Exercise: Using `rep`, how to generate 
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


Or, using the appropriate constructor and defining the number of columns and/or of rows:


```r
m <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
m <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, nrow = 3)
```


What happens if
- the number elements does not match?
- the number of cols/rows don't match the number of elements?

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


### Difference between `[` and `[[`

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

### Example

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
rm(age, consent, firstName, gender, secondName, weight)
```


By default, `character`s are converted to `factor`s when generating
`data.frame`s; use `stringsAsFactors = FALSE` to keep the `characters`
as is and convert to `factor`s explicitly.

Extract subsets of interest: patients that have given consent, male
patients, males that have given consent, those that are below the
average weight, order the `data.frame` by patient age, ...


```r
patients[patients$consent, ]
patients[patients$gender == "Male", ]
patients[patients$gender == "Male" & patients$consent, ]
```



```r
avgw <- mean(patients$weight)
patients[patients$weight < avgw, ]
```



```r
patients[order(patients$age), ]
```


## Inspecting data

- structure: `str()`
- `n` first elements or rows: `head()`
- `n` last elements or rows: `tail()`
- element names: `names()`, `colnames()`, `rownames()`
- size: `length()`, `dim()`, `nrow()`, `ncol()`
- `summary()`
- plotting!
- ...


## Exercise

Using the `weatherdata()` function from the `camweather` package,
download a weather data frame of your choice.
- What weather data is available? See `?weatherdata` and inspect the
  data frame's column names.
- What were the highest and lowest temperatures on that day? Hint: see
  `min`, `max` and/or `range` functions.
- The `Sun` and `Rain` values are cumulative from `Start`. What is the
  average rainfall per hours for that day? Hint: see `diff` for to
  calculate differences between successive values and `mean`.
- In what direction has the wind blown most on that day? Hint:
  `table`.

## A note of accessor speed

There are multiple ways to access elements in various `R` objects that
have different advantages. Using names, for instance, is very
convenient and avoids off-index errors, at the cost of accessor
timing. Below, we compare accessor timings (measured with
`system.time` and summed over 100 replications using the `replicate`
function).


```r
sum(replicate(100, system.time(mtcars["Volvo 142E", "carb"])["elapsed"]))
sum(replicate(100, system.time(mtcars[32, 11])["elapsed"]))
sum(replicate(100, system.time(mtcars[[11]][32])["elapsed"]))
sum(replicate(100, system.time(mtcars$carb[32])["elapsed"]))
```


See [here](https://gist.github.com/hadley/8150051) for a more detailed
comparison of the above example and
[here](https://gist.github.com/lgatto/8249301) for an example with
(named) lists.

## Collating matrices and data frames

It is easy to collate compatible matrices and data frames along their
columns or rows using `cbind` and `rbind`.


```r
m1 <- matrix(1:12, ncol = 3)
m2 <- matrix(1:9, ncol = 3)
m3 <- matrix(1:16, ncol = 4)
cbind(m1, m3)
rbind(m1, m2)
```


But


```r
cbind(m1, m2)
rbind(m1, m3)
```


## Special values: `NULL`, `NA`, `NaN`


```r
class(NA)
class(NaN)
class(NULL)
```

The `NULL` object:


```r
l[[2]] <- NULL
l
length(NULL)
c(1, NULL)
list(1, NULL)
```


But, `NA` can take different specific values for different atomic types:


```r
class(as.integer(NA))
class(as.numeric(NA))
class(as.character(NA))
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

## Environments

An `environment` is and un-ordered collection of symbols, or
associative arrays. They are implemented as hash tables.


```r
e <- new.env()
e
e$a <- 1
e$a
ls()  ## list content of global environment
ls(e)  ## list content of e
a <- 10  ## a different variable a
e$a
e[["a"]]
```


Values from specific environments can also be retrieved with `get` or
`mget` for multiple values or assigned with `assign`.

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

 
