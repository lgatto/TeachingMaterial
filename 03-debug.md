---
title: "Part III: Debugging"
author: "Laurent Gatto"
---

## Overview

- Debbugging: techniques and tools
- Condition handling: try/tryCatch
- Defensive programming
- Testing

## Defensive programming

Before we begin with debugging, let's look at ways to prevent bugs
(more at the end of this part). 

*Defensive programming*: 
- making the code work in a predicable manner
- writing code that fails in a well-defined manner
- if something *weird* happens, either properly deal with it, of fail
  quickly and loudly

The level of defensiveness will depend whether you write a function
for interactive of programmatic usage.

## Talk to your users

### Diagnostic messages


```r
message("This is a message for our dear users.")
```


```r
message("This is a message for our dear users. ",
	paste("Thank you for using our software",
              sw, "version", packageVersion(sw)))
```

Do not use `print` or `cat`:


```r
f1 <- function() {
    cat("I AM LOUD AND YOU CAN'T HELP IT.\n")
    ## do stuff
    invisible(TRUE)
}
f1()
```


```r
f2 <- function() {
    message("Sorry to interup, but...")
    ## do stuff
    invisible(TRUE)
}
f2()
suppressMessages(f2())
```

Of course, it is also possible to manually define verbosity. This
makes you write more code for a feature readily available. But still
better to use `message`.


```r
f3 <- function(verbose = TRUE) {
    if (verbose)
        message("I am being verbose because you let me.")
    ## do stuff
    invisible(TRUE)
}
f3()
f3(verbose = FALSE)
```

### Progress bars

- `utils::txtProgressBar` function


```r
n <- 10
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
    setTxtProgressBar(pb, i)
    Sys.sleep(0.5)
}
close(pb)
```

- [`progress`](https://github.com/gaborcsardi/progress) package


```r
library("progress")
pb <- progress_bar$new(total = n)
for (i in 1:n) {
    pb$tick()
    Sys.sleep(0.5)
}
```

Tip: do not over use progress bars. Ideally, a user should be
confident that everything is under control and progress is made while
waiting for a function to return. In my experience, a progress bar is
usefull when there is a specific and/or user-defined number of
iterations, such a *iterating over n files*, or *running a simulation
n times*.

**Question**: What about mixing progress bars and verbosity.

### Warning

> There is a problem with warnings. No one reads them. Pat Burns, in
> *R inferno*.


```r
warning("Do not ignore me. Somthing bad might have happened.")
warning("Do not ignore me. Somthing bad might be happening.", immediate. = TRUE)
```


```r
f <- function(...)
    warning("Attention, attention, ...!", ...)
f()
f(call. = FALSE)
```
Print warnings after they have been thrown.


```r
warnings()
last.warning
```

See also to `warn` option in `?options` .


```r
option("warn")
```

### Error


```r
stop("This is the end, my friend.")
```


```r
log(c(2, 1, 0, -1, 2)); print('end') ## warning 
xor(c(TRUE, FALSE));  print ('end')  ## error
```

Stop also has a `call.` parameter.


```r
geterrmessage()
```

### Logging

See for example the [`log4r`](https://github.com/johnmyleswhite/log4r)
package:


```r
## Import the log4r package.
library('log4r')

## Create a new logger object with create.logger().
logger <- create.logger()

## Set the logger's file output: currently only allows flat files.
logfile(logger) <- file.path('base.log')

## Set the current level of the logger.
level(logger) <- "INFO"

## Try logging messages at different priority levels.
debug(logger, 'A Debugging Message') ## Won't print anything
info(logger, 'An Info Message')
warn(logger, 'A Warning Message')
error(logger, 'An Error Message')
fatal(logger, 'A Fatal Error Message')
```

## KISS

Keep your functions simple and stupid (and short). 

## Failing fast and well

> Bounds errors are ugly, nasty things that should be stamped out
> whenever possible. One solution to this problem is to use the
> `assert` statement. The `assert` statement tells C++, "This can
> never happen, but if it does, abort the program in a nice way." One
> thing you find out as you gain programming experience is that things
> that can "never happen" happen with alarming frequency. So just to
> make sure that things work as they are supposed to, it’s a good idea
> to put lots of self checks in your program. -- Practical C++
> Programming, Steve Oualline, O'Reilly.


```r
if (!condition) stop(...)
```


```r
stopifnot(TRUE)
stopifnot(TRUE, FALSE)
```

For example to test input classes, lengths, ...


```r
f <- function(x) {
    stopifnot(is.numeric(x), length(x) == 1)
    invisible(TRUE)
}

f(1)
f("1")
f(1:2)
f(letters)
```

The [`assertthat`](https://github.com/hadley/assertthat) package:


```r
x <- "1"
library("assertthat")
stopifnot(is.numeric(x))
assert_that(is.numeric(x))
assert_that(length(x) == 2)
```

* `is.flag(x)`: is x `TRUE` or `FALSE`? (a boolean flag)
* `is.string(x)`: is x a length 1 character vector?
* `has_name(x, nm)`, `x %has_name% nm`: does `x` have component `nm`?
* `has_attr(x, attr)`, `x %has_attr% attr`: does `x` have attribute `attr`?
* `is.count(x)`: is x a single positive integer?
* `are_equal(x, y)`: are `x` and `y` equal?
* `not_empty(x)`: are all dimensions of `x` greater than 0?
* `noNA(x)`: is `x` free from missing values?
* `is.dir(path)`: is `path` a directory?
* `is.writeable(path)`/`is.readable(path)`: is `path` writeable/readable?
* `has_extension(path, extension)`: does `file` have given `extension`?
  
* `assert_that()` signal an error
* `see_if()` returns a logical value, with the error message as an attribute.
* `validate_that()` returns `TRUE` on success, otherwise returns the error as
  a string.

## Consistency and predictability

Reminder of the innteractive use vs programming examples: 
- `[` and `drop` 
- `sapply`, `lapply`, `vapply`

Remember also the concept of *tidy data*.

## Comparisons

### Floating point issues to be aware of

R FAQ [7.31](http://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f)?



```r
a <- sqrt(2)
a * a == 2
a * a - 2
```


```r
1L + 2L == 3L
1.0 + 2.0 == 3.0
0.1 + 0.2 == 0.3
```

### Floating point: how to compare

- `all.equal` compares R objects for *near equality*. Takes into
  account whether object attributes and names ought the taken into
  consideration (`check.attributes` and `check.names` parameters) and
  tolerance, which is machine dependent.


```r
all.equal(0.1 + 0.2, 0.3)
all.equal(0.1 + 0.2, 3.0)
isTRUE(all.equal(0.1 + 0.2, 3)) ## when you just want TRUE/FALSE
```

### Exact identity

`identical`: test objects for exact equality


```r
1 == NULL
all.equal(1, NULL)
identical(1, NULL)
identical(1, 1.)   ## TRUE in R (both are stored as doubles)
all.equal(1, 1L)
identical(1, 1L)   ## stored as different types
```

Appropriate within `if`, `while` condition statements. (not
`all.equal`, unless wrapped in `isTRUE`).

## Exercise

(From [adv-r](http://adv-r.had.co.nz/Exceptions-Debugging.html#defensive-programming).)

The `col_means` function computes the means of all numeric columns in
a data frame.


```r
col_means <- function(df) {
  numeric <- sapply(df, is.numeric)
  numeric_cols <- df[, numeric]
  data.frame(lapply(numeric_cols, mean))
}
```

Is it a robust function? What happens if there are unusual inputs.


```r
col_means(mtcars)
col_means(mtcars[, 0])
col_means(mtcars[0, ])
col_means(mtcars[, "mpg", drop = F])
col_means(1:10)
col_means(as.matrix(mtcars))
col_means(as.list(mtcars))

mtcars2 <- mtcars
mtcars2[-1] <- lapply(mtcars2[-1], as.character)
col_means(mtcars2)
```

## Debugging: techniques and tools

### Shit happens

> Funding your bug is a process of confirming the many things that you
> believe are true - until you find one which is not true. -- Norm Matloff

#### 1. Identify the bug (the difficult part)
- Something went wrong!
- Where in the code does it happen?
- Does it happen every time?
- What input triggered it?
- Report it (even if it is in your code - use github issues, for
  example).

**Tip**: Beware of your intuition. As a scientist, do what you are
used to: generate a hypotheses, *design an experiment* to test them,
and record the results.

#### 2. Fix it (the less difficult part)
- Correct the bug.
- Make sure that bug will not repeat itself!
- How can we be confident that we haven't introduced new bugs?

## Tools

- `print`/`cat`
- `traceback()`
- `browser()`
- `options(error = )`, `options(warn = )`
- `trace`
- RStudio

## Condition handling

- withCallingHandlers
- try/tryCatch


```r
safelog <- function(x) {
  tryCatch(log(x),
           error = function(e) "an error",
           warning = function(e) "a warning")
}
```

## Testing
