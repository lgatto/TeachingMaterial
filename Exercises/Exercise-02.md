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
##  [1]  3 13 15  1  6 11  2 17 12 14 16 19  7  8 10  9 18  4  5 20
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
##  [1] 0.19290 0.61971 0.48283 0.55939 0.47877 0.38178 0.01488 0.19948
##  [9] 0.54318 0.34219 0.50203 0.52519 0.60642 0.79278 0.93961 0.83412
## [17] 0.48865 0.09721 0.97350 0.32690
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
##  [1] "c" "m" "o" "a" "f" "k" "b" "q" "l" "n" "p" "s" "g" "h" "j" "i" "r"
## [18] "d" "e" "t"
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
##  [1] FALSE  TRUE  TRUE FALSE FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
## [12]  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE  TRUE
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



