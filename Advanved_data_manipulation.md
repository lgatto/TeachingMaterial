Bioinformatics for Big Omics Data: Advanced data manipulation
========================================================
width: 1440
height: 900
transition: none
font-family: 'Helvetica'
css: my_style.css
author: Raphael Gottardo, PhD
date: January 10, 2014

<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png" /></a><br /><tiny>This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US">Creative Commons Attribution-ShareAlike 3.0 Unported License</tiny></a>.

Motivation
==========

Let's first turn on the cache for increased performance.

```r
# Set some global knitr options
opts_chunk$set(cache=TRUE)
```



- R has pass-by-value semantics, which minimizes accidental side effects. However, this can become a major bottleneck when dealing with large datasets. 

- Fortunately, R provides some solution to this problems.


---------


```r
a <- rnorm(10^5)
gc()
```

```
         used (Mb) gc trigger (Mb) max used (Mb)
Ncells 314992 16.9     467875 25.0   407500 21.8
Vcells 563950  4.4    1031040  7.9   900303  6.9
```

```r
b <- a
gc()
```

```
         used (Mb) gc trigger (Mb) max used (Mb)
Ncells 315105 16.9     531268 28.4   407500 21.8
Vcells 564119  4.4    1031040  7.9   900303  6.9
```

```r
a[1] <- 0
gc()
```

```
         used (Mb) gc trigger (Mb) max used (Mb)
Ncells 315124 16.9     597831 32.0   407500 21.8
Vcells 664150  5.1    1162592  8.9   900303  6.9
```



Overview
========

Here we will review three R packages that can be used to provide efficient data manipulation:

- `data.table`: An package for efficient data storage and manipulation
- `RSQLite`: Database Interface R driver for SQLite
- `sqldf`: An R package for runing SQL statements on R data frames, optimized for convenience

<small>Thank to Kevin Ushey for the `data.table` notes</small>

What is data.table?
===================

`data.table` is an R package that extends. `R` `data.frame`s.

Under the hood, they are just `data.frame's, with some extra 'stuff' added on.
So, they're collections of equal-length vectors. Each vector can be of
different type.


```r
library(data.table)
dt <- data.table(x=1:3, y=c(4, 5, 6), z=letters[1:3])
dt
```

```
   x y z
1: 1 4 a
2: 2 5 b
3: 3 6 c
```

```r
class(dt)
```

```
[1] "data.table" "data.frame"
```


The extra functionality offered by `data.table` allows us to modify, reshape, 
and merge `data.table`s much quicker than `data.frame`s. **See that `data.table` inherits from `data.frame`!**

**Note:** The [development version (v1.8.11)](https://r-forge.r-project.org/scm/viewvc.php/pkg/NEWS?view=markup&root=datatable) of data.table includes a lot of new features including `melt` and `dcase` methods


Installing data.table
=====================

- stable CRAN release
    
    

```r
install.packages("data.table")
```

- latest bug-fixes + enhancements

```r
install.packages("data.table", repos="http://R-forge.R-project.org")
```



What's Different?
=================

Most of your interactions with `data.table`s will be through the subset (`[`)
operator, which behaves quite differently for `data.table`s. We'll examine
a few of the common cases.

Visit [this stackoverflow question](http://stackoverflow.com/questions/13618488/what-you-can-do-with-data-frame-that-you-cant-in-data-table) for a summary of the differences between `data.frame`s and `data.table`s.

Single element subsetting
=========================


```r
library(data.table)
DF <- data.frame(x=1:3, y=4:6, z=7:9)
DT <- data.table(x=1:3, y=4:6, z=7:9)
DF[c(2,3)]
```

```
  y z
1 4 7
2 5 8
3 6 9
```

```r
cat("\n")
```

```r
DT[c(2,3)]
```

```
   x y z
1: 2 5 8
2: 3 6 9
```


By default, single-element subsetting in `data.table`s refers to rows, rather
than columns.

Row subsetting
===============


```r
library(data.table)
DF <- data.frame(x=1:3, y=4:6, z=7:9)
DT <- data.table(x=1:3, y=4:6, z=7:9)
DF[c(2,3), ]
```

```
  x y z
2 2 5 8
3 3 6 9
```

```r
cat("\n")
```

```r
DT[c(2,3), ]
```

```
   x y z
1: 2 5 8
2: 3 6 9
```


Notice: row names are lost with `data.table`s. Otherwise, output is identical.

Column subsetting
=================


```r
library(data.table)
DF <- data.frame(x=1:3, y=4:6, z=7:9)
DT <- data.table(x=1:3, y=4:6, z=7:9)
DF[, c(2,3)]
```

```
  y z
1 4 7
2 5 8
3 6 9
```

```r
cat("\n")
```

```r
DT[, c(2,3)]
```

```
[1] 2 3
```


`DT[, c(2,3)]` just returns `c(2, 3)`. Why on earth is that?

The j Expression
================

The subset operator is really a function, and `data.table` modifies it to behave
differently.

Call the arguments we pass e.g. `DT[i, j]`, or `DT[i]`.

The second argument to `[` is called the `j expression`, so-called because it's
interpreted as an `R` expression. This is where most of the `data.table`
magic happens.

`j` is an expression evaluated within the frame of the `data.table`, so
it sees the column names of `DT`. Similarly for `i`.

First, let's remind ourselves what an `R` expression is.


Expressions
===========

An `expression` is a collection of statements, enclosed in a block generated by
braces `{}`.


```r
## an expression with two statements
{
  x <- 1
  y <- 2
}
## the last statement in an expression is returned
k <- { print(10); 5 }
```

```
[1] 10
```

```r
print(k)
```

```
[1] 5
```


The j Expression (suite)
================

So, `data.table` does something special with the `expression` that you pass as
`j`, the second argument, to the subsetting (`[`) operator.

The return type of the final statement in our expression determines the type
of operation that will be performed.

In general, the output should either be a `list` of symbols, or a statement
using `:=`.

We'll start by looking at the `list` of symbols as an output.

An example
==========

When we simply provide a `list` to the `j expression`, we generate a new
`data.table` as output, with operations as performed within the `list` call.


```r
library(data.table)
DT <- data.table(x=1:5, y=1:5)
DT[, list(mean_x = mean(x), sum_y = sum(y), sumsq=sum(x^2+y^2))]
```

```
   mean_x sum_y sumsq
1:      3    15   110
```


Notice how the symbols `x` and `y` are looked up within the `data.table` `DT`.
No more writing `DT$` everywhere!


Using :=
=========
Using the `:=` operator tells us we should assign columns by reference into
the `data.table` `DT`:


```r
library(data.table)
DT <- data.table(x=1:5)
DT[, y := x^2]
```

```
   x  y
1: 1  1
2: 2  4
3: 3  9
4: 4 16
5: 5 25
```

```r
print(DT)
```

```
   x  y
1: 1  1
2: 2  4
3: 3  9
4: 4 16
5: 5 25
```


Using := (suite)
================
By default, `data.table`s are not copied on a direct assignment `<-`:


```r
library(data.table)
DT <- data.table(x=1)
DT2 <- DT
DT[, y := 2]
```

```
   x y
1: 1 2
```

```r
DT2
```

```
   x y
1: 1 2
```


Notice that `DT2` has changed. This is something to be mindful of; if you want
to explicitly copy a `data.table` do so with `DT2 <- copy(DT)`.

A slightly more complicated example
===================================


```r
library(data.table)
DT <- data.table(x=1:5, y=6:10, z=11:15)
DT[, m := log2( (x+1) / (y+1) )]
```

```
   x  y  z       m
1: 1  6 11 -1.8074
2: 2  7 12 -1.4150
3: 3  8 13 -1.1699
4: 4  9 14 -1.0000
5: 5 10 15 -0.8745
```

```r
print(DT)
```

```
   x  y  z       m
1: 1  6 11 -1.8074
2: 2  7 12 -1.4150
3: 3  8 13 -1.1699
4: 4  9 14 -1.0000
5: 5 10 15 -0.8745
```



Using an expression in j
========================

Note that the right-hand side of a `:=` call can be an expression.


```r
library(data.table)
DT <- data.table(x=1:5, y=6:10, z=11:15)
DT[, m := { tmp <- (x + 1) / (y + 1); log2(tmp) }]
```

```
   x  y  z       m
1: 1  6 11 -1.8074
2: 2  7 12 -1.4150
3: 3  8 13 -1.1699
4: 4  9 14 -1.0000
5: 5 10 15 -0.8745
```

```r
print(DT)
```

```
   x  y  z       m
1: 1  6 11 -1.8074
2: 2  7 12 -1.4150
3: 3  8 13 -1.1699
4: 4  9 14 -1.0000
5: 5 10 15 -0.8745
```


Multiple returns from an expression in j
========================================

The left hand side of a `:=` call can also be a character vector of names,
for which the corresponding final statement in the `j expression` should be
a list of the same length.


```r
library(data.table)
DT <- data.table(x=1:5, y=6:10, z=11:15)
DT[, c('m', 'n') := { tmp <- (x + 1) / (y + 1); list( log2(tmp), log10(tmp) ) }]
```

```
   x  y  z       m       n
1: 1  6 11 -1.8074 -0.5441
2: 2  7 12 -1.4150 -0.4260
3: 3  8 13 -1.1699 -0.3522
4: 4  9 14 -1.0000 -0.3010
5: 5 10 15 -0.8745 -0.2632
```

```r
DT[, `:=`(a=x^2, b=y^2)]
```

```
   x  y  z       m       n  a   b
1: 1  6 11 -1.8074 -0.5441  1  36
2: 2  7 12 -1.4150 -0.4260  4  49
3: 3  8 13 -1.1699 -0.3522  9  64
4: 4  9 14 -1.0000 -0.3010 16  81
5: 5 10 15 -0.8745 -0.2632 25 100
```

```r
DT[, c("c","d"):=list(x^2, y^2)]
```

```
   x  y  z       m       n  a   b  c   d
1: 1  6 11 -1.8074 -0.5441  1  36  1  36
2: 2  7 12 -1.4150 -0.4260  4  49  4  49
3: 3  8 13 -1.1699 -0.3522  9  64  9  64
4: 4  9 14 -1.0000 -0.3010 16  81 16  81
5: 5 10 15 -0.8745 -0.2632 25 100 25 100
```


The j expression revisited 
===============

So, we typically call `j` the `j expression`, but really, it's either:

1. An expression, or

2. A call to the function `:=`, for which the first argument is a set of
names (vectors to update), and the second argument is an expression, with
the final statement typically being a list of results to assign within
the `data.table`.

As I said before, `a := b` is parsed by `R` as `":="(a, b)`, hence it
looking somewhat like an operator.


```r
quote(a := b)
```

```
`:=`(a, b)
```


Why does it matter?
===================

Whenever you sub-assign a `data.frame`, `R` is forced to copy the entire
`data.frame`.

That is, whenever you write `DF$x <- 1`, `DF["x"] <- 1`, `DF[["x"]] <- 1`...

... R will make a copy of `DF` before assignment.

This is done in order to ensure any other symbols pointing at the same object
do not get modified. This is a good thing for when we need to reason about
the code we write, since, in general, we expect `R` to operate without side 
effects.

Unfortunately, it is prohibitively slow for large objects, and hence why
`:=` can be very useful.

Why does it matter?
===================


```r
library(data.table); library(microbenchmark)
big_df <- data.frame(x=rnorm(1E6), y=sample(letters, 1E6, TRUE))
big_dt <- data.table(big_df)
microbenchmark( big_df$z <- 1, big_dt[, z := 1] )
```

```
Unit: milliseconds
                 expr   min      lq  median      uq     max neval
        big_df$z <- 1 19.69 105.380 107.212 109.875 467.243   100
 big_dt[, `:=`(z, 1)]  4.49   4.743   4.862   4.995   8.889   100
```


Once again, notice that `:=` is actually a function, and `z := 1` 
is parsed as `":="(z, 1)`.



Using by
========
We can also perform split-apply-combine like operations through the use of the
`by` argument:


```r
library(data.table)
DT <- data.table(x=1:5, y=6:10, gp=c('a', 'a', 'a', 'b', 'b'))
DT[, z := mean(x+y), by=gp]
```

```
   x  y gp  z
1: 1  6  a  9
2: 2  7  a  9
3: 3  8  a  9
4: 4  9  b 14
5: 5 10  b 14
```

```r
print(DT)
```

```
   x  y gp  z
1: 1  6  a  9
2: 2  7  a  9
3: 3  8  a  9
4: 4  9  b 14
5: 5 10  b 14
```


Notice that since `mean(x+y)` returns a scalar (numeric vector of length 1),
it is recycled to fill within each group.

Generating a new data.table
============================

What if, rather than modifying the current `data.table`, we wanted to generate
a new one?


```r
library(data.table)
DT <- data.table(x=1:5, y=6:10, gp=c('a', 'a', 'a', 'b', 'b'))
DT[, list(z=mean(x + y)), by=gp]
```

```
   gp  z
1:  a  9
2:  b 14
```


Notice that we retain one row for each unique group specified in the `by`
argument, and only the `by` variables along-side our `z` variable are returned.

The j Expression
================

- A `list`

... returns a new `data.table`, potentially subset over groups in your `by`.

- A `:=` Call

... modifies that `data.table` in place, hence saving memory. Output is
recycled if the `by` argument is used.

#### In general, our `j expression` is either:

2. an expression, with the final (or only) statement being a `list` of (named) 
  arguments,

3. a call to the `:=` function, with
  * the first argument being names, and
  * the second argument being an expression, for which the last statement is
  a list of the same length as the first argument.



Special Variables
=================

There are a number of special variables defined only within `j`, that allow us to
do some neat things...


```r
library(data.table)
data.table()[, ls(all=TRUE)]
```

```
[1] ".GRP"  ".I"    ".iSD"  ".N"    ".SD"   "print"
```


These variables allow us to infer a bit more about what's going on within the
`data.table` calls, and also allow us to write more complicated `j expression`s.


Special Variables
==================
## `.SD`

A `data.table` containing the subset of data for each group, excluding
columns used in `by`.

## `.BY`

A ` list` containing a length 1 vector for each item in `by`.

## `.N`

An integer, length 1, containing the number of rows in `.SD`.

## `.I`

A vector of indices, holding the row locations from which `.SD` was
pulled from the parent `DT`.

## `.GRP`

A counter telling you which group you're working with (1st, 2nd, 3rd...)

Example usage of .N - Counts
============================

Compute the counts, by group, using `data.table`...

```r
set.seed(123); library(data.table); library(microbenchmark)
DT <- data.table(x=sample(letters[1:3], 1E5, TRUE))
DT[, .N, by=x]
```

```
   x     N
1: a 33387
2: c 33201
3: b 33412
```

```r
table(DT$x)
```

```

    a     b     c 
33387 33412 33201 
```


Example usage of .N - Counts
============================



```r
library(data.table)
library(microbenchmark)
DT <- data.table(x=factor(sample(letters[1:3], 1E5, TRUE)))
microbenchmark( tbl=table(DT$x), DT=DT[, .N, by=x] )
```

```
Unit: milliseconds
 expr   min    lq median     uq   max neval
  tbl 6.897 8.980  9.620 13.224 21.41   100
   DT 2.840 3.319  3.448  4.147 14.79   100
```


Example usage of .SD - lapply-type calls
========================================


```r
library(data.table)
DT <- data.table(x=rnorm(10), y=rnorm(10), z=rnorm(10), id=letters[1:10])
DT[, sapply(.SD, mean), .SDcols=c('x', 'y', 'z')]
```

```
      x       y       z 
-0.2512  0.7089  0.4616 
```

```r
sapply(DT[,1:3, with=FALSE], mean)
```

```
      x       y       z 
-0.2512  0.7089  0.4616 
```


Example usage of .SD - lapply-type calls
========================================


```r
library(data.table); library(microbenchmark)
DT <- data.table(x=rnorm(1E5), y=rnorm(1E5), z=rnorm(1E5))
microbenchmark(
  DT=DT[, sapply(.SD, mean), .SDcols=c('x', 'y', 'z')],
  base=sapply(DT, mean)
)
```

```
Unit: milliseconds
 expr   min    lq median    uq    max neval
   DT 3.406 4.545  4.674 5.194 23.308   100
 base 2.386 2.443  2.492 2.720  4.215   100
```


`data.table` isn't always faster, unfortunately.

Keys
=====

`data.table`s can be keyed, allowing for faster indexing and subsetting. Keys
are also used for `join`s, as we'll see later.


```r
library(data.table)
DT <- data.table(x=c('a', 'a', 'b', 'c', 'a'), y=rnorm(5))
setkey(DT, x)
DT['a'] ## grabs rows corresponding to 'a'
```

```
   x      y
1: a 0.8760
2: a 0.6238
3: a 0.8211
```


Note that this does a `binary search` rather than a `vector scan`, which is
much faster!

Key performance
================


```r
library(data.table); library(microbenchmark)
DF <- data.frame(key=sample(letters, 1E6, TRUE), x=rnorm(1E6))
DT <- data.table(DF)
setkey(DT, key)
identical( DT['a']$x, DF[ DF$key == 'a', ]$x )
```

```
[1] TRUE
```

```r
microbenchmark( DT=DT['a'], DF=DF[ DF$key == 'a', ], times=5 )
```

```
Unit: milliseconds
 expr     min      lq  median      uq     max neval
   DT   1.805   1.998   2.189   2.226   2.582     5
   DF 143.524 144.564 151.459 151.743 159.640     5
```


Further reading: you can set multiple keys with `setkeyv` as well.

Joins
=====

`data.table` comes with many kinds of joins, implements through the 
`merge.data.table` function, and also through the `[` syntax as well. We'll
focus on using `merge`.


```r
library(data.table); library(microbenchmark)
DT1 <- data.table(x=c('a', 'a', 'b', 'dt1'), y=1:4)
DT2 <- data.table(x=c('a', 'b', 'dt2'), z=5:7)
setkey(DT1, x)
setkey(DT2, x)
merge(DT1, DT2)
```

```
   x y z
1: a 1 5
2: a 2 5
3: b 3 6
```


Overview of joins
=================

Here is a quick summary of SQL joins, applicable to `data.table` too.


<small>(Source: http://www.codeproject.com)</small>

![SQL joins](http://www.codeproject.com/KB/database/Visual_SQL_Joins/Visual_SQL_JOINS_orig.jpg)

A left join
===========


```r
library(data.table); library(microbenchmark)
DT1 <- data.table(x=c('a', 'a', 'b', 'dt1'), y=1:4)
DT2 <- data.table(x=c('a', 'b', 'dt2'), z=5:7)
setkey(DT1, x)
setkey(DT2, x)
merge(DT1, DT2, all.x=TRUE)
```

```
     x y  z
1:   a 1  5
2:   a 2  5
3:   b 3  6
4: dt1 4 NA
```


A right join
=============


```r
library(data.table); library(microbenchmark)
DT1 <- data.table(x=c('a', 'a', 'b', 'dt1'), y=1:4)
DT2 <- data.table(x=c('a', 'b', 'dt2'), z=5:7)
setkey(DT1, x)
setkey(DT2, x)
merge(DT1, DT2, all.y=TRUE)
```

```
     x  y z
1:   a  1 5
2:   a  2 5
3:   b  3 6
4: dt2 NA 7
```


An outer join
==============


```r
library(data.table); library(microbenchmark)
DT1 <- data.table(x=c('a', 'a', 'b', 'dt1'), y=1:4)
DT2 <- data.table(x=c('a', 'b', 'dt2'), z=5:7)
setkey(DT1, x)
setkey(DT2, x)
merge(DT1, DT2, all=TRUE) ## outer join
```

```
     x  y  z
1:   a  1  5
2:   a  2  5
3:   b  3  6
4: dt1  4 NA
5: dt2 NA  7
```


---

Speed example
=============


```r
library(data.table); library(microbenchmark)
DT1 <- data.table(
  x=do.call(paste, expand.grid(letters, letters, letters, letters)), 
  y=rnorm(26^4)
)
DT2 <- DT1[ sample(1:nrow(DT1), 1E5), ]
setnames(DT2, c('x', 'z'))
DF1 <- as.data.frame(DT1)
DF2 <- as.data.frame(DT2)
setkey(DT1, x); setkey(DT2, x)
microbenchmark( DT=merge(DT1, DT2), DF=merge.data.frame(DF1, DF2), replications=5)
```

```
Unit: nanoseconds
         expr       min        lq    median        uq       max neval
           DT 122458615 126900331 1.298e+08 1.520e+08 3.512e+08   100
           DF 858472417 931070468 9.856e+08 1.096e+09 1.227e+09   100
 replications         9        46 2.295e+02 3.000e+02 5.370e+02   100
```


Subset joins
============

We can also perform joins of two keyed `data.table`s using the `[` operator.
We perform right joins, so that e.g.

- `DT1[DT2]`

is a right join of `DT1` into `DT2`. These joins are typically a bit faster. Do
note that the order of columns post-merge can be different, though.

Subset Joins
============


```r
library(data.table); library(microbenchmark)
DT1 <- data.table(
  x=do.call(paste, expand.grid(letters, letters, letters, letters)), 
  y=rnorm(26^4)
)
DT2 <- DT1[ sample(1:nrow(DT1), 1E5), ]
setnames(DT2, c('x', 'z'))
setkey(DT1, x); setkey(DT2, x)
tmp1 <- DT2[DT1]
setcolorder(tmp1, c('x', 'y', 'z'))
tmp2 <- merge(DT1, DT2, all.x=TRUE)
setcolorder(tmp2, c('x', 'y', 'z'))
identical(tmp1, tmp2)
```

```
[1] TRUE
```


Subset Joins can be Faster
==========================


```r
library(data.table); library(microbenchmark)
DT1 <- data.table(
  x=do.call(paste, expand.grid(letters, letters, letters, letters)), 
  y=rnorm(26^4)
)
DT2 <- DT1[ sample(1:nrow(DT1), 1E5), ]
setnames(DT2, c('x', 'z'))
setkey(DT1, x); setkey(DT2, x)
microbenchmark( bracket=DT1[DT2], merge=merge(DT1, DT2, all.y=TRUE), times=5 )
```

```
Unit: milliseconds
    expr    min    lq median     uq    max neval
 bracket  39.58  40.9  42.34  42.38  46.75     5
   merge 220.40 221.0 248.67 287.39 308.82     5
```


data.table and SQL
===================

We can understand the usage of `[` as SQL statements. 
From [FAQ 2.16](http://datatable.r-forge.r-project.org/datatable-faq.pdf):

data.table Argument | SQL Statement
 ---|---
i | WHERE
j | SELECT
:= | UPDATE
by | GROUP BY
i | ORDER BY (in compound syntax)
i | HAVING (in compound syntax)

Compound syntax refers to multiple subsetting calls, and generally isn't
necessary until you really feel like a `data.table` expert:

    DT[where,select|update,group by][having][order by][ ]...[ ]

data.table and SQL - Joins
===================

Here is a quick summary table of joins in `data.table`.


SQL | data.table
 ---|---
LEFT JOIN | x[y]
RIGHT JOIN | y[x]
INNER JOIN | x[y, nomatch=0]
OUTER JOIN | merge(x,y)

data.table and SQL
==================

It's worth noting that I really mean it when I say that `data.table` is like
an in-memory data base. It will even perform some basic query optimization!


```r
library(data.table)
options(datatable.verbose=TRUE)
DT <- data.table(x=1:5, y=1:5, z=1:5, a=c('a', 'a', 'b', 'b', 'c'))
DT[, lapply(.SD, mean), by=a]
```

```
Finding groups (bysameorder=FALSE) ... done in 0secs. bysameorder=TRUE and o__ is length 0
Optimized j from 'lapply(.SD, mean)' to 'list(.External(Cfastmean, x, FALSE), .External(Cfastmean, y, FALSE), .External(Cfastmean, z, FALSE))'
Starting dogroups ... done dogroups in 0 secs
```

```
   a   x   y   z
1: a 1.5 1.5 1.5
2: b 3.5 3.5 3.5
3: c 5.0 5.0 5.0
```



Some thoughts
============

The primary use of `data.table` is for the fast split-apply-combine functionality.


```r
library(data.table); library(microbenchmark)
DT <- data.table(gp1=sample(letters, 1E6, TRUE), gp2=sample(LETTERS, 1E6, TRUE), y=rnorm(1E6))
microbenchmark( times=5,
  DT=DT[, mean(y), by=list(gp1, gp2)],
  DF=with(DT, tapply(y, paste(gp1, gp2), mean))
)
```

```
Finding groups (bysameorder=FALSE) ... done in 0.072secs. bysameorder=FALSE and o__ is length 1000000
Detected that j uses these columns: y 
Optimized j from 'mean(y)' to '.External(Cfastmean, y, FALSE)'
Starting dogroups ... done dogroups in 0.062 secs
Finding groups (bysameorder=FALSE) ... done in 0.099secs. bysameorder=FALSE and o__ is length 1000000
Detected that j uses these columns: y 
Optimized j from 'mean(y)' to '.External(Cfastmean, y, FALSE)'
Starting dogroups ... done dogroups in 0.075 secs
Finding groups (bysameorder=FALSE) ... done in 0.081secs. bysameorder=FALSE and o__ is length 1000000
Detected that j uses these columns: y 
Optimized j from 'mean(y)' to '.External(Cfastmean, y, FALSE)'
Starting dogroups ... done dogroups in 0.073 secs
Finding groups (bysameorder=FALSE) ... done in 0.118secs. bysameorder=FALSE and o__ is length 1000000
Detected that j uses these columns: y 
Optimized j from 'mean(y)' to '.External(Cfastmean, y, FALSE)'
Starting dogroups ... done dogroups in 0.055 secs
Finding groups (bysameorder=FALSE) ... done in 0.104secs. bysameorder=FALSE and o__ is length 1000000
Detected that j uses these columns: y 
Optimized j from 'mean(y)' to '.External(Cfastmean, y, FALSE)'
Starting dogroups ... done dogroups in 0.046 secs
```

```
Unit: milliseconds
 expr   min    lq median    uq   max neval
   DT 155.4 180.6  181.4 188.5 359.0     5
   DF 583.8 665.1  755.0 901.8 957.2     5
```


Other interesting convenience functions
===============

- `like`


```r
DT = data.table(Name=c("Mary","George","Martha"), Salary=c(2,3,4))
# Use regular expressions
DT[Name %like% "^Mar"]
```

```
     Name Salary
1:   Mary      2
2: Martha      4
```


- `set*` functions
`set`, `setattr`, `setnames`, `setcolorder`, `setkey`, `setkeyv`


```r
setcolorder(DT, c("Salary", "Name"))
DT
```

```
   Salary   Name
1:      2   Mary
2:      3 George
3:      4 Martha
```



- `DT[, (myvar):=NULL]` remove a column


```r
DT[,Name:=NULL]
```

```
   Salary
1:      2
2:      3
3:      4
```



Bonuses: fread
===============
`data.table` also comes with `fread`, a file reader much, much better than
`read.table` or `read.csv`:


```r
library(data.table); library(microbenchmark)
big_df <- data.frame(x=rnorm(1E6), y=rnorm(1E6))
file <- tempfile()
write.table(big_df, file=file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
microbenchmark( fread=fread(file), r.t=read.table(file, header=TRUE, sep="\t"), times=1 )
```

```
Unit: milliseconds
  expr     min      lq  median      uq     max neval
 fread   289.6   289.6   289.6   289.6   289.6     1
   r.t 10120.0 10120.0 10120.0 10120.0 10120.0     1
```

```r
unlink(file)
```


Bonuses: rbindlist
==================

Use this function to `rbind` a list of `data.frame`s, `data.table`s or `list`s:


```r
library(data.table); library(microbenchmark)
dfs <- replicate(100, data.frame(x=rnorm(1E4), y=rnorm(1E4)), simplify=FALSE)
all.equal( rbindlist(dfs), data.table(do.call(rbind, dfs)) )
```

```
[1] TRUE
```

```r
microbenchmark( DT=rbindlist(dfs), DF=do.call(rbind, dfs), times=5 )
```

```
Unit: milliseconds
 expr     min      lq   median       uq      max neval
   DT   6.992   7.098    7.235    7.241    7.313     5
   DF 994.886 999.373 1000.934 1060.417 1209.355     5
```


Learning More
=============

- Read some of the `[data.table]` tagged questions on 
[StackOverflow](http://stackoverflow.com/questions/tagged/data.table)

- Read through the [data.table FAQ](http://datatable.r-forge.r-project.org/datatable-faq.pdf),
which is surprisingly well-written and comprehensive.

- Experiment!

Databases and the Structured Query Language (SQL) 
=================================================

- A database is an organized collection of datasets (tables).
- A database management system (DBMS) is a software system designed to allow the definition, creation, querying, update, and administration of databases.
- Well-known DBMSs include MySQL, PostgreSQL, SQLite, Microsoft SQL Server, Oracle, etc.
- Relational DBMSs (RDBMs) store data in a set of related tables
- Most RDBMs use some form of the Structured Query Language (SQL)

The Structured Query Language (SQL) 
===================================

Although SQL is an ANSI (American National Standards Institute) standard, there are different flavors of the SQL language.


The data in RDBMS is stored in database objects called tables.

A table is a collection of related data entries and it consists of columns and rows.

Here we will use SQLite, which is a self contained relational database management system. In contrast to other database management systems, SQLite is not a separate process that is accessed from the client application (e.g. MySQL, PostgreSQL).

Using RSQLite
===============================

Here we will make use of the [Bioconductor](http://www.bioconductor.org) project to load and use an SQLite database.


```r
# You only need to run this once
source("http://bioconductor.org/biocLite.R")
biocLite(c("org.Hs.eg.db"))
```



```r
# Now we can use the org.Hs.eg.db to load a database
library(org.Hs.eg.db)
# Create a connection
Hs_con <- org.Hs.eg_dbconn()
# List tables
head(dbListTables(Hs_con))
```

```
[1] "accessions"            "alias"                 "chrlengths"           
[4] "chromosome_locations"  "chromosomes"           "cytogenetic_locations"
```

```r
# Or using an SQLite command (NOTE: This is specific to SQLite)
head(dbGetQuery(Hs_con, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;"))
```

```
                   name
1            accessions
2                 alias
3            chrlengths
4  chromosome_locations
5           chromosomes
6 cytogenetic_locations
```


Using RSQLite (suite)
===============================


```r
# What columns are available?
dbListFields(Hs_con, "gene_info")
```

```
[1] "_id"       "gene_name" "symbol"   
```

```r
dbListFields(Hs_con, "alias")
```

```
[1] "_id"          "alias_symbol"
```

```r
# Or using SQLite
dbGetQuery(Hs_con, "PRAGMA table_info('gene_info');")
```

```
  cid      name         type notnull dflt_value pk
1   0       _id      INTEGER       1       <NA>  0
2   1 gene_name VARCHAR(255)       1       <NA>  0
3   2    symbol  VARCHAR(80)       1       <NA>  0
```



```r
gc()
```

```
          used (Mb) gc trigger (Mb) max used (Mb)
Ncells 1134968 60.7    1476915 78.9  1368491 73.1
Vcells 1153087  8.8    1925843 14.7  1623758 12.4
```

```r
alias <- dbGetQuery(Hs_con, "SELECT * FROM alias;")
gc()
```

```
          used (Mb) gc trigger (Mb) max used (Mb)
Ncells 1234655 66.0    1710298 91.4  1368491 73.1
Vcells 1504690 11.5    2287241 17.5  1883389 14.4
```

```r
gene_info <- dbGetQuery(Hs_con, "SELECT * FROM gene_info;")
chromosomes <- dbGetQuery(Hs_con, "SELECT * FROM chromosomes;")
```


-------------


```r
CD154_df <- dbGetQuery(Hs_con, "SELECT * FROM alias a JOIN gene_info g ON g._id = a._id WHERE a.alias_symbol LIKE 'CD154';")
gc()
```

```
          used (Mb) gc trigger  (Mb) max used (Mb)
Ncells 1278724 68.3    1967602 105.1  1368491 73.1
Vcells 1946078 14.9    2685683  20.5  2283449 17.5
```

```r
CD40LG_alias_df <- dbGetQuery(Hs_con, "SELECT * FROM alias a JOIN gene_info g ON g._id = a._id WHERE g.symbol LIKE 'CD40LG';")
gc()
```

```
          used (Mb) gc trigger  (Mb) max used (Mb)
Ncells 1278756 68.3    1967602 105.1  1368491 73.1
Vcells 1946189 14.9    2899967  22.2  2283449 17.5
```


Some SQL Commands
============

### SELECT

The SELECT is used to query the database and retrieve selected data that match the specific criteria that you specify:

SELECT column1 [, column2, ...] 
FROM tablename 
WHERE condition 

### ORDER BY

ORDER BY clause can order column name in either ascending (ASC) or descending (DESC) order.

### JOIN

There are times when we need to collate data from two or more tables. 
As with data.tables we can use LEFT/RIGHT/INNER JOINS

### GROUP BY

The GROUP BY was added to SQL so that aggregate functions could return a result grouped by column values. 

SELECT col_name, function (col_name) FROM table_name GROUP BY col_name 

A "GROUP BY" example
==================


```r
dbGetQuery(Hs_con, "SELECT c.chromosome, COUNT(g.gene_name) AS count FROM chromosomes c JOIN gene_info g ON g._id = c._id WHERE c.chromosome IN (1,2,3,4,'X') GROUP BY c.chromosome ORDER BY count;")
```

```
  chromosome count
1          4  1846
2          X  2134
3          3  2457
4          2  3021
5          1  4329
```


Some more SQL commands
=======================

Some other SQL statements that might be of used to you:

### CREATE TABLE

The CREATE TABLE statement is used to create a new table. 

### DELETE

The DELETE command can be used to remove a record(s) from a table. 

### DROP

To remove an entire table from the database use the DROP command.

### CREATE VIEW

A view is a virtual table that is a result of SQL SELECT statement. A view contains fields from one or more real tables in the database. This virtual table can then be queried as if it were a real table. 


Creating your own SQLite database in R
===================================


```r
db <- dbConnect(SQLite(), dbname="./Data/SDY61/SDY61.sqlite")
dbWriteTable(conn = db, name = "hai", value = "./Data/SDY61/hai_result.txt", row.names = FALSE, header = TRUE, sep="\t")
```

```
[1] FALSE
```

```r
dbWriteTable(conn = db, name = "cohort", value = "./Data/SDY61/arm_or_cohort.txt", row.names = FALSE, header = TRUE, sep="\t")
```

```
[1] FALSE
```



Creating your own SQLite database in R (suite)
=============================================


```r
dbListFields(db, "hai")
```

```
 [1] "RESULT_ID"                 "ARM_ACCESSION"            
 [3] "BIOSAMPLE_ACCESSION"       "EXPSAMPLE_ACCESSION"      
 [5] "EXPERIMENT_ACCESSION"      "STUDY_ACCESSION"          
 [7] "STUDY_TIME_COLLECTED"      "STUDY_TIME_COLLECTED_UNIT"
 [9] "SUBJECT_ACCESSION"         "VALUE_REPORTED"           
[11] "VIRUS_STRAIN"              "WORKSPACE_ID"             
```

```r
dbListFields(db, "cohort")
```

```
[1] "ARM_ACCESSION"             "DESCRIPTION"              
[3] "NAME"                      "POPULATION_SELECTION_RULE"
[5] "SORT_ORDER"                "STUDY_ACCESSION"          
[7] "TYPE"                      "WORKSPACE_ID"             
```

```r

dbGetQuery(db, "SELECT STUDY_TIME_COLLECTED, cohort.DESCRIPTION, MAX(VALUE_REPORTED) AS max_value FROM hai JOIN cohort ON hai.ARM_ACCESSION = cohort.ARM_ACCESSION WHERE cohort.DESCRIPTION LIKE '%TIV%' GROUP BY BIOSAMPLE_ACCESSION;")
```

```
   STUDY_TIME_COLLECTED                      DESCRIPTION max_value
1                     0 Healthy adults given TIV vaccine        80
2                     0 Healthy adults given TIV vaccine        40
3                     0 Healthy adults given TIV vaccine       160
4                     0 Healthy adults given TIV vaccine       160
5                     0 Healthy adults given TIV vaccine         5
6                     0 Healthy adults given TIV vaccine        40
7                     0 Healthy adults given TIV vaccine        40
8                     0 Healthy adults given TIV vaccine        40
9                     0 Healthy adults given TIV vaccine      2560
10                   28 Healthy adults given TIV vaccine       640
11                   28 Healthy adults given TIV vaccine       640
12                   28 Healthy adults given TIV vaccine       160
13                   28 Healthy adults given TIV vaccine       640
14                   28 Healthy adults given TIV vaccine      2560
15                   28 Healthy adults given TIV vaccine       160
16                   28 Healthy adults given TIV vaccine      2560
17                   28 Healthy adults given TIV vaccine       640
18                   28 Healthy adults given TIV vaccine      2560
19                    0 Healthy adults given TIV vaccine         5
20                    0 Healthy adults given TIV vaccine       160
21                    0 Healthy adults given TIV vaccine       640
22                    0 Healthy adults given TIV vaccine        20
23                    0 Healthy adults given TIV vaccine        20
24                    0 Healthy adults given TIV vaccine       320
25                    0 Healthy adults given TIV vaccine         5
26                    0 Healthy adults given TIV vaccine        20
27                    0 Healthy adults given TIV vaccine       320
28                    0 Healthy adults given TIV vaccine         5
29                    0 Healthy adults given TIV vaccine        40
30                    0 Healthy adults given TIV vaccine       160
31                    0 Healthy adults given TIV vaccine       320
32                    0 Healthy adults given TIV vaccine       640
33                    0 Healthy adults given TIV vaccine       640
34                    0 Healthy adults given TIV vaccine        80
35                    0 Healthy adults given TIV vaccine        40
36                    0 Healthy adults given TIV vaccine         5
37                    0 Healthy adults given TIV vaccine        20
38                    0 Healthy adults given TIV vaccine         5
39                    0 Healthy adults given TIV vaccine         5
40                    0 Healthy adults given TIV vaccine         5
41                    0 Healthy adults given TIV vaccine         5
42                    0 Healthy adults given TIV vaccine         5
43                    0 Healthy adults given TIV vaccine        20
44                    0 Healthy adults given TIV vaccine        80
45                    0 Healthy adults given TIV vaccine        10
46                    0 Healthy adults given TIV vaccine       320
47                   28 Healthy adults given TIV vaccine       160
48                   28 Healthy adults given TIV vaccine       160
49                   28 Healthy adults given TIV vaccine       640
50                   28 Healthy adults given TIV vaccine        20
51                   28 Healthy adults given TIV vaccine       640
52                   28 Healthy adults given TIV vaccine       320
53                   28 Healthy adults given TIV vaccine        20
54                   28 Healthy adults given TIV vaccine        20
55                   28 Healthy adults given TIV vaccine       320
56                   28 Healthy adults given TIV vaccine         5
57                   28 Healthy adults given TIV vaccine       320
58                   28 Healthy adults given TIV vaccine       320
59                   28 Healthy adults given TIV vaccine       640
60                   28 Healthy adults given TIV vaccine       640
61                   28 Healthy adults given TIV vaccine       160
62                   28 Healthy adults given TIV vaccine       320
63                   28 Healthy adults given TIV vaccine       160
64                   28 Healthy adults given TIV vaccine        40
65                   28 Healthy adults given TIV vaccine        40
66                   28 Healthy adults given TIV vaccine       320
67                   28 Healthy adults given TIV vaccine        80
68                   28 Healthy adults given TIV vaccine        40
69                   28 Healthy adults given TIV vaccine       160
70                   28 Healthy adults given TIV vaccine       160
71                   28 Healthy adults given TIV vaccine        80
72                   28 Healthy adults given TIV vaccine        80
73                   28 Healthy adults given TIV vaccine        80
74                   28 Healthy adults given TIV vaccine        80
```


The sqldf package
=================

Sometimes it can be convenient to use SQL statements on dataframes. This is exactly what the sqldf package does.



```r
library(sqldf)
data(iris)
sqldf("select * from iris limit 5")
sqldf("select count(*) from iris")
sqldf("select Species, count(*) from iris group by Species")
```


The `sqldf` package can even provide increased speed over pure R operations. 

Summary
=======

- R base `data.frames` are convenient but often not adapted to large dataset manipulation (e.g. genomics). 

- Thankfully, there are good alternatives. My recommendtion is:
    - Use `data.table` for your day-to-day operations
    - When you have many tables and a complex schema, use `sqlite`.
    
**Note:** There many other R packages for "big data" such the `bigmemory` suite, `biglm`, `ff`, `RNetcdf`, `rhdf5`, etc. 
