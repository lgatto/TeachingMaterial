# Why test?

* Nobody writes buggy code... sure!

* How do you know your code is right unless you test?

* Can act as specification/examples of what the programme should do.

# Inverting some data can have serious consequences ...

http://www.sciencemag.org/content/314/5807/1856.full.pdf

2001 Science paper, and two more Science papers retracted after
inversion along x axis noted.

# Testing needs to be:

* reliable
* fast
* easy to run
* easy to summarise

# stop() vs warning()

* A warning is softer than an error; if a warning is generated
your program will still continue, whereas an error will stop the
program.


```r
log(c(2, 1, 0, -1, 2)); print('end')  # warning 
xor( c(TRUE, FALSE));  print ('end')  # error
```

* If you try to isolate warnings, you can change warnings to
  errors: `options(warn=2)`

* Add warnings and errors to your code using `warning()`,
  `stop()`


```r
x <- -1:1
y <- -5:-1

if ( any(x<0) ) warning("some elements of x are negative")
```

```
## Warning: some elements of x are negative
```

```r
if ( all(y<0) ) stop("all elements of x are negative")
```

```
## Error in eval(expr, envir, enclos): all elements of x are negative
```

```r
stopifnot(all.equal(pi, 3.1415927), 2 < 2, all(1:10 < 12), "a" < "b")
```

```
## Error: 2 < 2 is not TRUE
```

# Catching errors and warnings


```r
safelog <- function(x) {
  tryCatch(log(x),
           error = function(e) "an error",
           warning = function(e) "a warning")
}

safelog(3)
safelog(-5)
safelog("string")
```

# Floating point issues to be aware of

R FAQ 7.31? http://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f



```r
a <- sqrt(2)
a * a == 2
```

```
## [1] FALSE
```

```r
a * a - 2
```

```
## [1] 4.440892e-16
```


```r
1.0 + 2.0 == 3.0
0.1 + 0.2 == 0.3
```

# Floating point: how to compare


```r
all.equal( 0.1 + 0.2, 0.3)
```

```
## [1] TRUE
```

```r
all.equal( 0.1 + 0.2, 3.0)
```

```
## [1] "Mean relative difference: 9"
```

```r
isTRUE(all.equal(0.1 + 0.2, 3))           # when you just want TRUE/FALSE
```

```
## [1] FALSE
```

# Compile this document.


```r
library("knitr")
knit('testing.Rmd')
knit2html('testing.Rmd')
browseURL('testing.html')
```
