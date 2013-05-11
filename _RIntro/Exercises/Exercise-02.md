## A vector of `integer`s

We are going to use the `sample` function to create such a vector. 
As its name implied, the function samples among a set it is provided as input. 
It is as such important to provide a set in integers. As explained in 
`?sample`, if given a single numeric input `n`, `sample` will sample `n` values 
from `1:n` without replacement (`replacement = FALSE`), resulting in a permutation of `1:n`.

Using `n` equal `20` and making sure that `1:20` produces a sequence of integers


```r
class(1:20)
```

```
## [1] "integer"
```


we can now produce our vector of integers `i`:


```r
i <- sample(20)
i
```

```
##  [1]  2  7  5 17 14 10 15  1 12 13  6  9 18 11 16  3 19 20  4  8
```

```r
typeof(i)
```

```
## [1] "integer"
```


## A vector of `double`s

The `runif(n)` function samples `n` values from a uniform distribution ranging (by default) from 0 to 1. 


```r
d <- runif(20)
d
```

```
##  [1] 0.43514 0.28883 0.87758 0.08802 0.32376 0.15818 0.76988 0.51394
##  [9] 0.93103 0.14682 0.62064 0.34712 0.40829 0.44959 0.80837 0.78791
## [17] 0.16037 0.28768 0.79473 0.48247
```

```r
typeof(d)
```

```
## [1] "double"
```


We could also have used `rnorm` to sample values from a normal distribution. 


```r
typeof(rnorm(20))
```

```
## [1] "double"
```

## A vector of `character`s

Instead of typing the sequence of characters ourselves, we will use the 
built-in variable `letters`, which stores the 26 letters from the alphabet. 
We will select the letters at positions 1 to 20 in the order defined by `i`.


```r
head(letters)
```

```
## [1] "a" "b" "c" "d" "e" "f"
```

```r
s <- letters[i]
s
```

```
##  [1] "b" "g" "e" "q" "n" "j" "o" "a" "l" "m" "f" "i" "r" "k" "p" "c" "s"
## [18] "t" "d" "h"
```

```r
typeof(s)
```

```
## [1] "character"
```

```r
length(s)
```

```
## [1] 20
```


Other built-in variables of interest are `LETTERS` (as `letters`, but capital letters), 
`month.name` (month names), pi, ... see `?letters`.

## A vector of `logical`s

Although `TRUE` is a perfectly valid vector of logicals (of length 1), we are going to use 
a comparison operator on the vector `i` to create a vector of logicals 
(see `?'=='` or `help("<")` for details). 
Each element of `i` will be used in turn and compared to the second operand (`10` below), 
effectively resulting in a vector of logicals of the same length than `i`.


```r
l <- i >= 10
## same as
l <- (i >= 10)
l
```

```
##  [1] FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE
## [12] FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE FALSE
```

```r
typeof(l)
```

```
## [1] "logical"
```

```r
length(i)
```

```
## [1] 20
```

```r
length(l)
```

```
## [1] 20
```



