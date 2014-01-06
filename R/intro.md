Introduction
===




# Pre-requisites

## Basic tools

**Variables** are used to store values: they bind a value to a name. The
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
`install.packages(, repos = NULL)`. However, this will force you to
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


or, if you have downloaded the package source
([`camweather_0.1.2.tar.gz`](http://proteome.sysbiol.cam.ac.uk/lgatto/files/camweather_0.1.2.tar.gz))
or windows binary
([`camweather_0.1.2.zip`](http://proteome.sysbiol.cam.ac.uk/lgatto/files/camweather_0.1.2.zip))


```r
install.packages("camweather_0.1.2.zip", repos = NULL)  ## on windows
install.packages("camweather_0.1.2.tar.gz", repos = NULL)  ## elsewhere
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
```

```
## [1] "y" "z"
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
help(package = "devtools")
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
```

```
## R Under development (unstable) (2013-10-16 r64064)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] knitr_1.5
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.1 formatR_0.10   stringr_0.6.2  tools_3.1.0
```

```r
version
```

```
##                _                                                 
## platform       x86_64-unknown-linux-gnu                          
## arch           x86_64                                            
## os             linux-gnu                                         
## system         x86_64, linux-gnu                                 
## status         Under development (unstable)                      
## major          3                                                 
## minor          1.0                                               
## year           2013                                              
## month          10                                                
## day            16                                                
## svn rev        64064                                             
## language       R                                                 
## version.string R Under development (unstable) (2013-10-16 r64064)
## nickname       Unsuffered Consequences
```


# `R` is

<!-- TODO: write intro paragraph -->

## A dynamic language


```r
a <- 1L
typeof(a)
```

```
## [1] "integer"
```

```r

b <- 1
typeof(b)
```

```
## [1] "double"
```

```r

typeof(a + b)
```

```
## [1] "double"
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

## Object-oriented

With multiple frameworks: S3, S4 and S4 reference classes. See
[R object oriented programming](https://github.com/lgatto/roo),
[Short S4 tutorial](https://github.com/lgatto/S4-tutorial) and a
[A (Not So) Short Introduction to S4](http://cran.r-project.org/doc/contrib/Genolini-S4tutorialV0-5en.pdf)
for more details.

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

All these can be initialised with their respective constructor
functions:


```r
vector()
matrix()
array()
list()
data.frame()
```


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


### Factors


```r
gender <- factor(c("male", "female", "male"))
gender
class(gender)
typeof(gender)
```


These specific vector types can also be initialised specifically:


```r
character()
logical()
numeric()
double()
factor()
```


## Vectorised operations


```r
(x <- c(1, 2, 3))
```

```
## [1] 1 2 3
```

```r
(y <- c(3, 2, 1))
```

```
## [1] 3 2 1
```

```r
x + y
```

```
## [1] 4 4 4
```

```r
x^2
```

```
## [1] 1 4 9
```

```r
sqrt(x)
```

```
## [1] 1.000 1.414 1.732
```



```r
x <- c("a", "b", "c")
paste(x, 1, sep = ".")
```

```
## [1] "a.1" "b.1" "c.1"
```


### Recycling


```r
x <- c(1, 2, 3, 4)
y <- c(1, 2)
x + y
```

```
## [1] 2 4 4 6
```

```r
z <- c(1, 2, 3)
x + z
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## [1] 2 4 6 5
```



```r
x <- c(1, 2, 3, 4)
x[2:3] <- 10
x
```

```
## [1]  1 10 10  4
```



## Generating vectors

- `seq` and **argument matching by position and name**
- `:`
- **Exercise**: Using `rep`, how to generate 
  - 3 repetitions of 1 to 10: 1, 2, ..., 9, 10, 1, 2, ..., 9, 10, 1, 2, ..., 9, 10
  - repeating numbers 1 to 10 each 3 times: 1, 1, 1, 2, 2, 2, ..., 9, 9, 9, 10, 10, 10
- `rnorm`, `runif`, ... to draw values from specific distributions.
- `sample` to create permutations of vectors. 

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
## or provideDimnames(m, base = list(letters, LETTERS))
m
m["A", "c"]
```


Dimensions can also have their own names:


```r
M <- matrix(c(4, 8, 5, 6, 4, 2, 1, 5, 7), nrow=3)
dimnames(M) <- list(year =
                    c(2005, 2006, 2007),
                    "mode of transport" =
                    c("plane", "bus", "boat"))
M
```


## Simplification/dropping of dimensions


```r
m[-1, ]  ## becomes a vector
m[-1, , drop = FALSE]  ## remains a matrix
m[, -(1:2)]  ## becomes a vector
m[, -(1:2), drop = FALSE]  ## remains a matrix
```

In `R`, there is no such thing as a columns/row vector:


```r
x <- 1:6
t(x)  ## becomes a matrix
dim(t(x))
dim(t(t(x)))
```


## Array

Like a matrix with > 2 dimensions. 


```r
x <- 1:30
dim(x) <- c(5, 3, 2)
## or array(1:30, dim = c(5, 3, 2))
x
```


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
`names` (same as `colnames`) and `rownames` or indices can be
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


### Example


```r
chrsms13 <- weatherdata("2013-12-25")
chrsms12 <- weatherdata("2012-12-25")
chrsms11 <- weatherdata("2011-12-25")
chrsms <- rbind(chrsms11, chrsms12, chrsms13)
nrow(chrsms11) + nrow(chrsms12) + nrow(chrsms13)
nrow(chrsms)
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

 
