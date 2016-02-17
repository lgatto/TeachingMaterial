---
title: "Part I: Introduction"
author: "Laurent Gatto"
---

## References

- [Advanced R](http://adv-r.had.co.nz/), Hadley Wickham.
- [The R Inferno](http://www.burns-stat.com/documents/books/the-r-inferno/), Patrick Burns.
- [An Introduction to the Interactive Debugging Tools in R](http://www.biostat.jhsph.edu/~rpeng/docs/R-debug-tools.pdf), Roger D. Peng.
- [R Programming for Bioinformatics](http://master.bioconductor.org/help/publications/books/r-programming-for-bioinformatics/), Robert Gentleman.

## Overview

- Coding style(s)
- Interactive use and programming
- Environments
- Computing on the language
- Tidy data


## Introduction

> Computers are cheap, and thinking hurts. -- Use Ligges

Simplicity, readability and consistency are a long way towards
robust code.

## Coding style(s)

Why?

> Good coding style is like using correct punctuation. You can manage
> without it, but it sure makes things easier to read.
-- Hadley Wickham

for **consistency** and **readability**.

## Which one?

- [Bioconductor](http://master.bioconductor.org/developers/how-to/coding-style/)
- [Hadley Wickham](http://r-pkgs.had.co.nz/style.html)
- [Google](http://google.github.io/styleguide/Rguide.xml)
- ...

## Examples

- Place spaces around all infix operators (`=`, `+`, `-`, `<-`, etc., but *not* `:`)
  and after a comma (`x[i, j]`).
- Spaces before `(` and after `)`.
- Use `<-` rather than `=`.
- Limit your code to 80 characters per line
- Indentation: do not use tabs, use 2 (HW)/4 (Bioc) spaces
- Function names: use verbs
- Variable names: camelCaps (Bioc)/ `_` (HW) (but not a `.`)
- Prefix non-exported functions with a ‘.’ (Bioc).
- Class names: start with a capital
- Comments: `# ` or `## ` (from emacs)

## [`formatR`](https://cran.rstudio.com/web/packages/formatR/index.html)


```r
library("formatR")
tidy_eval(text = c("a=1+1;a  # print the value", "matrix ( rnorm(10),5)"),
          arrow = TRUE)
```

```
## a <- 1 + 1
## a  # print the value
## ## [1] 2
## 
## matrix(rnorm(10), 5)
## ##            [,1]       [,2]
## ## [1,] -0.1366365  1.7361156
## ## [2,] -0.2361091 -0.5469096
## ## [3,] -0.7019853 -0.6514337
## ## [4,] -0.4159157  0.4958135
## ## [5,] -0.7996917 -0.9635818
```

## [`BiocCheck`](http://bioconductor.org/packages/devel/bioc/html/BiocCheck.html)

```
* Checking function lengths................
  The longest function is 677 lines long
  The longest 5 functions are:
* Checking formatting of DESCRIPTION, NAMESPACE, man pages, R source,
  and vignette source...
    * CONSIDER: Shortening lines; 616 lines (11%) are > 80 characters
      long.
    * CONSIDER: Replacing tabs with 4 spaces; 3295 lines (60%) contain
      tabs.
    * CONSIDER: Indenting lines with a multiple of 4 spaces; 162 lines
      (2%) are not.
```

## Style changes over time

![Style changes over time](./figs/style.png)

## Interactive use vs programming: `drop`


```r
head(cars)
head(cars[, 1])
head(cars[, 1, drop = FALSE])
```

## Interactive use vs programming: `sapply/lapply`

```
df1 <- data.frame(x = 1:3, y = LETTERS[1:3])
sapply(df1, class)
df2 <- data.frame(x = 1:3, y = Sys.time() + 1:3)
sapply(df2, class)
```
## Ineractive use vs programming

Moving from using R to programming R is *abstraction*, *automation*,
*generalisation*.

## Semantics

- *pass-by-value* copy-on-modify
- *pass-by-reference*: environments, S4 Reference Classes

## Environments

### Motivation

- Data structure that enables *scoping* (see later).
- Have reference semantics
- Useful data structure on their own

### Definition (1)

An environment associates, or *binds*, names to values in memory.
Variables in an environment are hence called *bindings*.

## Creating and populate environments


```r
e <- new.env()
e$a <- 1
e$b <- LETTERS[1:5]
e$c <- TRUE
e$d <- mean
```


```r
e$a <- e$b
e$a <- LETTERS[1:5]
```

- Objects in environments have unique names
- Objects in different environments can of course have identical names
- Objects in an environment have no order
- Environments have parents

## Definition (2)

An environment is composed of a *frame* that contains the name-object
bindings and a parent (enclosing) environment.

## Relationship between environments

Every environment has a parent (enclosing) environment


```r
e <- new.env()
parent.env(e)
```
Current environment


```r
environment()
```

```
## <environment: R_GlobalEnv>
```

Noteworthy environments


```r
globalenv()
emptyenv()
baseenv()
```

All parent of `R_GlobalEnv`:


```r
search()
```

```
## [1] ".GlobalEnv"        "package:formatR"   "package:stats"    
## [4] "package:graphics"  "package:grDevices" "package:utils"    
## [7] "package:datasets"  "Autoloads"         "package:base"
```

```r
as.environment("package:stats")
```

```
## <environment: package:stats>
## attr(,"name")
## [1] "package:stats"
## attr(,"path")
## [1] "/usr/local/lib64/R/library/stats"
```

Listing objects in an environment


```r
ls() ## default is R_GlobalEnv
```

```
## character(0)
```

```r
ls(envir = e)
```

```
## Error in ls(envir = e): object 'e' not found
```

```r
ls(pos = 1)
```

```
## character(0)
```


```r
search()
```

```
## [1] ".GlobalEnv"        "package:formatR"   "package:stats"    
## [4] "package:graphics"  "package:grDevices" "package:utils"    
## [7] "package:datasets"  "Autoloads"         "package:base"
```

Note: Every time a package is loaded with `library`, it is inserted in
the search path after the `R_GlobalEnv`.

## Accessors and setters

- In addition to `$`, one can also use `[[`, `get` and `assign`.
- To check if a name exists in an environmet (or in any or its
  parents), one can use `exists`.
- Compare two environments with `identical` (not `==`).

## Exercise

- Draw a few environments with variables and ask to reproduce in R.

## Where is a symbol defined?

`pryr::where()`

## Lexical scoping

[Lexical comes from *lexical analysis* in computer science, which is the conversion of characters (code) into a sequence of meaningful (for the computer) tokens.]

**Definition**: Rules that define how R looks up values for a given name/symbol.

- Objects in environments have unique names
- Objects in different environments can of course have identical names.
- If a name is not found in the current environment, it is looked up
  in the parent (enclosing) from. 
- If it is not found in the parent (enclosing) frame, it is looked up
  in the parent's parent frame, and so on...


```r
search()
mean <- function(x) cat("The mean is", sum(x)/length(x), "\n")
mean(1:10)
base::mean(1:10)
rm(mean)
mean(1:10)
```

## Exercises

Start by mentally running the code chunks below - what do the functions return? 

After testing new code chunks, don't forget to clean up your
workspace, to avoid unexpected results.


```r
f <- function() {
    x <- 1
    y <- 2
    c(x, y)
}
f()
```


```r
x <- 2
g <- function(){
    y <- 1
    c(x, y)
}
g()
```


```r
x <- 1
h <- function() {
    y <- 2
    i <- function() {
        z <- 3
        c(x, y, z)
    }
    i()
}
h()
```


```r
j <- function(x) {
    y <- 2
    function(){
        c(x, y)
    }
}
k <- j(1)
k()
```


```r
j <- function() {
    if (!exists("a")) {
        a <- 1
    } else {
        a <- a + 1
    }
    print(a)
}
j() ## First call
j() ## Secong call
```


```r
f <- function() x
x <- 1
f() 
x <- 2
f()
```


```rm
f <- function(x) {
    f <- function(x) {
        f <- function(x) {
            x^2
        }
        f(x) + 1
    }
    f(x) * 2
}
f(10)
```

## Scoping

(See also later functions)

*Lexical scoping*: default behaviour, current environment, then
traversing *enclosing/parent environments*.

*Dynamic scoping*: looking up variables in the *calling environment*,
used in non-standard evaluation.

## Assignments

- `<-` assigns/creates in the current environment

- `<<-` (deep assignment) never creates a variable in the current
  environment, but modifies an existing variable in the current or
  first enclosing environment where that name is defined. 
  
- If `<<-` does not find the name, it will create the variable in the global environment.

## Using environments

Most environments are created when creating and calling
functions. They are also used in packages:

- Used in packages: *package* and *namespace* environments

There are several reasons to create then manually.

- Reference semantics
- Avoiding copies
- Package state
- As a hashmap for fast name lookup

## Reference semantics


```r
modify <- function(x) {
	x$a <- 2
	invisible(TRUE)
}
```


```r
x_l <- list(a = 1)
modify(x_l)
x_l$a
```


```r
x_e <- new.env()
x_e$a <- 1
modify(x_e)
x_e$a
```

Tip: when setting up environments, it is advised to set to parent
(enclosing) environment to be `emptyenv()`, to avoid accidentally
inheriting objects from somewhere else on the search path.


```r
e <- new.env(parent.env = empty.env())
```

```
## Error in new.env(parent.env = empty.env()): unused argument (parent.env = empty.env())
```

### Exercise

What is going to happen when we access `"x"` in the four cases below?


```r
x <- 1
e1 <- new.env()
get("x", envir = e1)
```


```r
get("x", envir = e1, inherits = FALSE)
```


```r
e2 <- new.env(parent = emptyenv())
get("x", envir = e2)
```


```r
get("x", envir = e1, inherits = FALSE)
```

## Avoiding copies

Since environments have reference semantics, they are not copied.
When passing an environment as function argument (directly, or as part
of a more complex data structure), it is **not** copied: all its
values are accessible within the function and can be persistently
modified.


```r
e <- new.env()
e$x <- 1
f <- function(myenv) myenv$x <- 2
f(e)
e$x
```

This is used in the `eSet` class family to store the expression data.


```r
library("Biobase")
getClass("eSet")
getClass("AssayData")
new("ExpressionSet")
```

## Preserving state in packages

Explicit envirionments are also useful to preserve state or define
constants-like variables in a package. One can then set getters and
setters for users to access the variables within that private
envionment.

#### Use case

Colour management in [`pRoloc`](https://github.com/lgatto/pRoloc/blob/master/R/environment.R):


```r
.pRolocEnv <- new.env(parent=emptyenv(), hash=TRUE)

stockcol <- c("#E41A1C", "#377EB8", "#238B45", "#FF7F00", "#FFD700", "#333333",
              "#00CED1", "#A65628", "#F781BF", "#984EA3", "#9ACD32", "#B0C4DE",
              "#00008A", "#8B795E", "#FDAE6B", "#66C2A5", "#276419", "#CD8C95",
              "#6A51A3", "#EEAD0E", "#0000FF", "#9ACD32", "#CD6090", "#CD5B45",
              "#8E0152", "#808000", "#67000D", "#3F007D", "#6BAED6", "#FC9272")

assign("stockcol", stockcol, envir = .pRolocEnv)

getStockcol <- function() get("stockcol", envir = .pRolocEnv)

setStockcol <- function(cols) {
    if (is.null(cols)) {
        assign("stockcol", stockcol, envir = .pRolocEnv)
    } else { 
		assign("stockcol", cols, envir = .pRolocEnv)
	}
}
```

and in plotting functions:


```r
...
if (missing(col))
  col <- getStockcol()
...
```

Hadley's tip: Invisibly returning the old value from 


```r
setStockcol <- function(cols) {
	prevcols <- getStockcol()
    if (is.null(cols)) {
        assign("stockcol", stockcol, envir = .pRolocEnv)
    } else { 
		assign("stockcol", cols, envir = .pRolocEnv)
	}
	invisible(prevcols)
}
```

## Computing on the language

- character containing variables names
- use arguments to name things
- loading data

## Tidy data

(interactive use)

> Hadley Wickham, Tidy Data, Vol. 59, Issue 10, Sep 2014, Journal of
> Statistical Software. http://www.jstatsoft.org/v59/i10.

Tidy datasets are easy to manipulate, model and visualize, and have a
specific structure: each variable is a column, each observation is a
row, and each type of observational unit is a table.


## Tidy tools

Tidy data also makes it easier to develop tidy tools for data
analysis, tools that both input and output tidy datasets.

- `dply::select` select columns
- `dlpy::filter` select rows
- `dplyr:mutate` create new columns
- `dpplyr:group_by` split-apply-combine
- `dlpyr:summarise` collapse each group into a single-row summary of
  that group
- `magrittr::%>%` piping


## Examples


```r
library("dplyr")
surveys <- read.csv("http://datacarpentry.github.io/dc_zurich/data/portal_data_joined.csv")
head(surveys)

surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight)
  
surveys %>%
  mutate(weight_kg = weight / 1000) %>%
  filter(!is.na(weight)) %>%
  head

surveys %>%
  group_by(sex) %>%
  tally()
  
surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))
  
surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE),
            min_weight = min(weight, na.rm = TRUE)) %>%
  filter(!is.nan(mean_weight))
```

## Application to other data structures

> Hadley Wickham (@hadleywickham) tweeted at 8:45 pm on Fri, Feb 12,
> 2016: @mark_scheuerell @drob **the importance of tidy data is not the
> specific form, but the consistency**
> (https://twitter.com/hadleywickham/status/698246671629549568?s=09)

- Well-formatted and well-documented `S4` class 
- `S4` as input -function-> `S4` as output

![MSnSet schematics](https://raw.githubusercontent.com/lgatto/pRoloc/master/vignettes/Figures/msnset.png)
