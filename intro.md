# Introduction


## Writing R functions

A function is made of
- a name
- some inputs (formal parameters)
- a single output (return value)
- a body
- an environment, the map of the location of the functions variable

```{r, eval=FALSE}
f <- function(x) {
    y <- x + 1
    return(x * y)
}
```

## The `apply` functions

A functional approach to iteration:

```
apply(X, MARGIN, FUN, ...)
```

- `MARGIN` = 1 for row, 2 for cols.
- `FUN` = function to apply
- `...` = extra args to function.
- `simplify` =  should the result be simplified if possible.

To iterate over vectors and lists, there's `sapply` and `lapply`,
where the margin is implicit.

```
lapply(X, FUN, ...) 
sapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) 
```
