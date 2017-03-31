---
title: "Part I: Introduction"
author: "Laurent Gatto"
---

## Overview

- Coding style(s)
- Interactive use and programming
- Environments

## Introduction

> Computers are cheap, and thinking hurts. -- Uwe Ligges

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
- Spaces before `(` and after `)`; not for function.
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
tidy_source(text = "a=1+1;a  # print the value
                    matrix ( rnorm(10),5)",
            arrow = TRUE)
```

```
## a <- 1 + 1
## a  # print the value
## matrix(rnorm(10), 5)
```

## [`BiocCheck`](http://bioconductor.org/packages/devel/bioc/html/BiocCheck.html)

```
$ R CMD BiocCheck package_1.0.0.tgz
```

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


## Ineractive use vs programming

Moving from using R to programming R is *abstraction*, *automation*,
*generalisation*.

## Interactive use vs programming: `drop`


```r
head(cars)
head(cars[, 1])
head(cars[, 1, drop = FALSE])
```

## Copying semantics: environments

Unlike other object types in R, environments have pass-by-reference semantics.
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

