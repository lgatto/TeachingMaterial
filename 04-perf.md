---
title: "Part IV: Performance"
author: "Laurent Gatto"
---

# Overview

- Benchmarking
- Profiling
- Optimisation
- Memory
- Rcpp
- Big data

# A word of caution

Knuth, Donald. /Structured Programming with `go to` Statements/, ACM
Journal Computing Surveys, Vol 6, No. 4, Dec. 1974. p.268.

> We should forget about small efficiencies, say about 97% of the
> time: premature optimization is the root of all evil. Yet we should
> not pass up our opportunities in that critical 3%. A good programmer
> will not be lulled into complacency by such reasoning, he will be
> wise to look carefully at the critical code; but only after that
> code has been identified


Robert Gentleman, in R Programming for Bioinformatics, 2008,
about R's built-in C interfaces:

> Since R is not compiled, in some situations its performance can be
> substantially improved by writing code in a compiled language. There
> are also reasons not to write code in other languages, and in
> particular we caution against premature optimization, prototyping in
> R is often cost effective. And in our experience very few routines
> need to be implemented in other languages for efficiency
> reasons. Another substantial reason not to use an implementation in
> some other language is increased complexity. The use of another
> language almost always results in higher maintenance costs and less
> stability. In addition, any extensions or enhancements of the code
> will require someone that is proficient in both R and the other
> language.

(`Rcpp` does make some of the above caution statements slightly less
critical.)

# R performance

R is not a fast language, but it is most of the time *fast enough* for
what we want to do, in particular with respect to interactive data
analysis. In such cases a slow but expressive and flexible language is
way better than a fast but less expressive and flexible
alternative. It is also relatively easy to avoid bad R programming
idioms that make code too slow.

## Timing, benchmarking

Let's compare two implementation of the square root calculation:
`sqrt(x)` and `x ^ 0.5`.


```r
x <- runif(100)
system.time(sqrt(x))
system.time(x^0.5)
```

Does this work? 


```r
x <- runif(1e5)
system.time(sqrt(x))
system.time(x^0.5)
```

We want to repeat timings multiple times:


```r
summary(replicate(10, system.time(x^0.5)[["elapsed"]]))
```

A better approach for such cases is the
[`microbenchmark` package](https://cran.rstudio.com/web/packages/microbenchmark/index.html)),
which is ideal to accurately benchmark small pieces of code, in
particular sub-millisecond (nanoseconds) executions (see units below). 

Each expression is run 100 times (controlled by the `times`
argument). In addition, the execution order is randomised and summary
timings are reported.


```r
x <- runif(100)
library(microbenchmark)

microbenchmark(sqrt(x),
               x ^ 0.5)
```

**Question**: where does this difference come from?

<!-- ```{r, eval=FALSE} -->
<!-- sqrt ## 1 argument -->
<!-- `^`  ## 2 arguments -->
<!-- as.list(body(function(x) x^0.5))   ## 2 symbols -->
<!-- as.list(body(function(x) sqrt(x))) ## 3 symbols -->
<!-- ``` -->

## Extreme dynamism

## Name lookup

## Lazy evaluation

## Extracting a single value from a data frame

## Improving R's performance

R implementation:

- [GNU R](http://www.r-project.org/)
- [pqR](http://www.pqr-project.org)
- [Renjin](http://www.renjin.org)
- [FastR](https://github.com/allr/fastr)
- [Riposte](https://github.com/jtalbot/riposte)
- [CXXR](http://www.cs.kent.ac.uk/project/cxxr/)

Example of *deferred evaluation*.

# Profiling


### Reminder

> We should forget about small efficiencies, say about 97% of the
> time: premature optimization is the root of all evil. Yet we should
> not pass up our opportunities in that critical 3%. A good programmer
> will not be lulled into complacency by such reasoning, he will be
> wise to look carefully at the critical code; but only after that
> code has been identified

Knuth, Donald. /Structured Programming with `go to` Statements/, ACM
Journal Computing Surveys, Vol 6, No. 4, Dec. 1974. p.268.

### Before optimising

1. Find the major bottleneck: code *profiling*.
2. Try to eliminate it.
3. Repeat until *fast enough*: ideally, define fast enough in advance.


### Make sure the code remains correct


```r
x <- runif(100)
all.equal(sqrt(x), x ^ 0.5)
```

```
## [1] TRUE
```
and unit tests.

### Are implementations really equivalent?


```r
library("sequences")
gccount
gccountr <- function(x) table(strsplit(x, "")[[1]])
gccountr2 <- function(x) tabulate(factor(strsplit(x, "")[[1]]))
```

Checking that our different implementations give the same results:


```r
s <- paste(sample(c("A", "C", "G", "T"),
                  100, replace = TRUE),
           collapse = "")

gccount(s)
gccountr(s)
gccountr2(s)
```

But are they really the same? Are we really comparing the same
functionalities?

Is it worth it?


```r
library("microbenchmark")
microbenchmark(gccount(s),
                     gccountr(s),
                     gccountr2(s),
                     times = 1e4, 
					 unit = "eps")
```


```r
library("ggplot2")
mb <- microbenchmark(gccount(s),
                     gccountr(s),
                     gccountr2(s))
print(mb)
microbenchmark:::autoplot.microbenchmark(mb)
```

## Profiling tools

- `Rprof` and `summaryRprof` functions: records timings at fixed
  intervals (default `interval` is 0.02 seconds)
- [`proftools`](https://cran.rstudio.com/web/packages/proftools/index.html)
  package:
- [`profr`](https://cran.rstudio.com/web/packages/profr/index.html)
  package:
- [`lineprof`](https://cran.rstudio.com/web/packages/lineprof/index.html)
  package: each *line of code* is profiled. This is less precise (than
  `Rprof`) but easier to interprete. Code must be sourced with
  `source()`.

# Optimisation

- look for existing solutions
- Do as little as possible: `gccountr` vs `gccountr2` example above;
  avoid names when not needed (`unlist`); simpler data.structures; use
  any assumptions about the data, at the cost of
  generalisation. Trade-off fast vs dangerous,
  flexibility/functionality vs performance.
- vectorisation: see  *R-vectorisation*
- parallelisation: see *R-parallel*
- avoid copies
- byte-code compilation: see section in *R-bmark-prof-optim*
- Re-implementing the code in C/C++ (see below)
- Case study: t-test in *Advanced R* - compare with `genefilter::rowttests`.

# Memory

- `tracemem`
- `object.size`
- `lineprof`

# Rcpp

See [here](https://github.com/lgatto/2016-02-25-adv-programming-EMBL/blob/master/rc.md).

# Big data

- [CRAN High-Performance and Parallel Computing task view](http://cran.r-project.org/web/views/HighPerformanceComputing.html).
- Storing data in database or databases-like structures: `RMySQL`,
      `RdbiPgSQL`, \ldots, `RSQLite`, `qldf`, `data.table` (the
      `data.table::fread`, when `read.table` is slow, also `scan`),
      `dplyr`, ... packages
- The `ff` package by Adler et al. offers file-based access to data
  sets that are too large to be loaded into memory, along with a
  number of higher-level functions
- The `bigmemory` package by Kane and Emerson permits storing large
  objects such as matrices in memory (as well as via files) and uses
  `external pointer` objects to refer to them
- `netCDF` data files: `ncdf` and `RNetCDF` packages
- `hdf5` format: `rhdf5` package
- `mmap` memory-mapped files/devices I/O
- hadoop and R
- See http://r-pbd.org/ and the
  [pbdDemo](http://cran.r-project.org/web/packages/pbdDEMO/)
  package/vignette.
- [Bioconductor in the cloud](http://bioconductor.org/help/bioconductor-cloud-ami/)
- [Bioconductor docker containers](http://bioconductor.org/help/docker/)
- ...

