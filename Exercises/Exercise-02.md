## A vector of `integer`s

We are going to use the `sample` function to create such a vector. 
As its name implied, the function samples among a set it is provided as imput. 
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
##  [1]  1  9 19  3 12 13  8  5  6  7 17 16 18 20 11  4  2 14 15 10
```

```r
typeof(i)
```

```
## [1] "integer"
```


## A vector of `double`s

The `runing` sample a certain number of values from a uniform distribution ranging (by default)
from 0 to 1. 


```r
d <- runif(20)
d
```

```
##  [1] 0.613792 0.590519 0.148436 0.337720 0.005709 0.772950 0.650761
##  [8] 0.072411 0.802297 0.131191 0.567193 0.440654 0.544043 0.666007
## [15] 0.811204 0.678165 0.336886 0.715581 0.342947 0.941654
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
##  [1] "a" "i" "s" "c" "l" "m" "h" "e" "f" "g" "q" "p" "r" "t" "k" "d" "b"
## [18] "n" "o" "j"
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
a logical operator on the vector `i` to create a vector of logicals. Each element of `i` 
will be used to the comparison in turn, effectively resulting in a vector of logicals of the 
same length than `i`.


```r
l <- i >= 10
l
```

```
##  [1] FALSE FALSE  TRUE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE
## [12]  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE
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



