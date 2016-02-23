---
title: "Part II: Functional programming"
author: "Laurent Gatto"
---

## Content

- Side effects
- Functions
- Functions environments
- More about functions
- Closures
- Higher order functions


## Writing functions

> To understand compuations in R, two slogans are helpful:
> - Everything that exists is an object.
> - Everything that happens is a function call.
> John Chambers

A function is made of
- a name
- some inputs (formal paramters)
- a single output (return value)
- a body
- an environment, the map of the location of the functions variable


```r
f <- function(x) {
    y <- x + 1
    return(x * y)
}
```

And these can be accessed and modified indivdually


```r
body(f)
args(f)
environment(f)

body(f) <- quote({
    y <- x * y
    return(x + y)
})
```

## Important

![Messy code hides bugs](./figs/funs.png) 

- Functions are a means of **abstraction**. A concept/computation is
  encapsulated/isolated from the rest with a function.
- Functions should **do one thing**, and do it well (compute, or plot,
  or save, ... not all in one go).
- **Side effects**: your functions should not have any (unless, of
  course, that is the main point of that function - plotting, write to
  disk, ...). Functions shouldn't make any changes in any
  environment. The only return their output.
- **Do not use global variables**. Everything the function needs is
  being passed as an argument. Function must be **self-contained**.
- Function streamline code and process

From the `R Inferno`:

Make your functions as simple as possible. Simple has many advantages:

- Simple functions are likely to be human efficient: they will be easy
  to understand and to modify.
- Simple functions are likely to be computer efficient.
- Simple functions are less likely to be buggy, and bugs will be
  easier to fix.
- (Perhaps ironically) simple functions may be more general—thinking
  about the heart of the matter often broadens the application.


Functions can be

1. Correct.
2. An error occurs that is clearly identified.
3. An obscure error occurs.
4. An incorrect value is returned.

We like category 1. Category 2 is the right behavior if the inputs do
not make sense, but not if the inputs are sensible. Category 3 is an
unpleasant place for your users, and possibly for you if the users
have access to you. Category 4 is by far the worst place to be - the
user has no reason to believe that anything is wrong. Steer clear of
category 4.


Finally, functions are

- Easier to debug (part III)
- Easier to profile (part IV)
- Easier to parallelise (part IV)

Functions are an central part of robust R programming.

## Lexical scoping

- If a name is not found in a functions environment, it is looked up
  in the parent (enclosing) from.
- If it is not found in the parent (enclosing) frame, it is looked up
  in the parent's parent frame, and so on...

*Lexical scoping*: default behaviour, current environment, then
traversing *enclosing/parent environments*.


```r
f <- function(x) x + y

f(1)

environment(f)
y <- 2
f(1)
```


```r
e <- new.env()
environment(f) <- e

f(1)
e$y <- 10
f(1)
```

This is of course bad practice, we don't want to rely on global variables.


```r
codetools::findGlobals(f)
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


```r
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

## More about functions

- Argument matching by position or by names
- Calling a function with a list of arguments


```r
args <- list(x = 1:10, trim = 0.3)
do.call(mean, args)
```

- Default arguments


```r
f <- function(x = 1, y = 2) x * y
f <- function(x = 1, y = x + 2) x * y
```

- Missing arguments


```r
f <- function(x = 1, y) {
	c(misssing(x), missing(y))
}
```

- Passing non-matched parameters `...` to an inner function


```r
plot2 <- function(...) {
    message("Verbose plotting...")
    plot(...)
}

f <- function(...) list(...)
```

- Return values: last statement, explicit `return`, make output
  `invivisble`
  

```r
f1 <- function() 1
f2 <- function() return(1)
f3 <- function() return(invisible(1))
```

- Explicit triggers before exiting. Useful to restore global state
  (plotting parameters, cleaning temporary files, ...)


```r
f1 <- function(x) {
    on.exit(print("!"))
    x + 1
}

f2 <- function(x) {
    on.exit(print("!"))
    stop("Error")
}
```


```r
f3 <- function() {
    on.exit(print("1"))
    on.exit(print("2"))
    invisible(TRUE)
}


f4 <- function() {
    on.exit(print("1"))
    on.exit(print("2"), add = TRUE)
    invisible(TRUE)
}
```


## More about scoping and environments

### Scoping

*Lexical scoping*: default behaviour, current environment, then
traversing *enclosing/parent environments*.


*Dynamic scoping*: looking up variables in the *calling environment*,
used in non-standard evaluation.



## More about functions

## Closures

## Higher order functions
