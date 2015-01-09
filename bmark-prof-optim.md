Benchmarking, profiling and optimisation
===============



## Warning

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

## Timing


```r
m <- matrix(runif(1e4), nrow=1000)
system.time(apply(m, 1, sum))
```

```
##    user  system elapsed 
##   0.001   0.000   0.002
```


```r
replicate(5, system.time(apply(m, 1, sum))[[1]])
```

```
## [1] 0.002 0.001 0.002 0.001 0.001
```


## Benchmarking


```r
library("sequences")
gccount
```

```
## function (inseq) 
## {
##     .Call("gccount", inseq, PACKAGE = "sequences")
## }
## <environment: namespace:sequences>
```

```r
gccountr <- function(x) table(strsplit(x, "")[[1]])
gccountr2 <- function(x) tabulate(factor(strsplit(x, "")[[1]]))
```

Checking that our different implementations give the same results:


```r
s <- paste(sample(c("A", "C", "G", "T"),
                  100, replace = TRUE),
           collapse = "")

gccount(s)
```

```
## [1] 23 33 25 19
```

```r
gccountr(s)
```

```
## 
##  A  C  G  T 
## 23 33 25 19
```

```r
gccountr2(s)
```

```
## [1] 23 33 25 19
```

But are they really the same? Are we really comparing the same
functionalities?


```r
library("microbenchmark")

mb <- microbenchmark(gccount(s),
                     gccountr(s),
                     gccountr2(s),
                     times = 1e4)
print(mb)
```

```
## Unit: microseconds
##          expr     min       lq       mean   median       uq       max
##    gccount(s)   2.101   2.8300   4.641646   5.1680   5.3850  1546.457
##   gccountr(s) 134.382 140.0335 148.821964 141.9350 145.5930  3243.746
##  gccountr2(s)  69.992  74.1940  82.524462  75.6785  77.6735 39558.145
##  neval cld
##  10000 a  
##  10000   c
##  10000  b
```


```r
library("ggplot2")
microbenchmark:::autoplot.microbenchmark(mb)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 


```r
library("rbenchmark")

benchmark(replications = 1e4,
          gccount(s),
          gccountr(s),
          gccountr2(s),
          columns=c('test', 'elapsed', 'replications'))          
```

```
##           test elapsed replications
## 3 gccountr2(s)   0.809        10000
## 2  gccountr(s)   1.577        10000
## 1   gccount(s)   0.051        10000
```

## Profiling 


```r
Rprof("rprof")
replicate(5, res <- apply(m, 1, mean, trim=.3))
Rprof(NULL)
summaryRprof("rprof")
```

```
$by.self
         self.time self.pct total.time total.pct
"apply"       0.04    33.33       0.12    100.00
"FUN"         0.04    33.33       0.12    100.00
"any"         0.02    16.67       0.02     16.67
"unique"      0.02    16.67       0.02     16.67

$by.total
               total.time total.pct self.time self.pct
"apply"              0.12    100.00      0.04    33.33
"FUN"                0.12    100.00      0.04    33.33
"lapply"             0.12    100.00      0.00     0.00
"replicate"          0.12    100.00      0.00     0.00
"sapply"             0.12    100.00      0.00     0.00
"mean.default"       0.04     33.33      0.00     0.00
"sort.int"           0.04     33.33      0.00     0.00
"any"                0.02     16.67      0.02    16.67
"unique"             0.02     16.67      0.02    16.67

$sample.interval
[1] 0.02

$sampling.time
[1] 0.12
```

### The `lineprof` package

The example below is from the
[`lineprof`](https://github.com/hadley/lineprof) github repo.


```r
library(lineprof)
source(find_ex("read-delim.r"))
wine <- find_ex("wine.csv")

x <- lineprof(read_delim(wine, sep = ","), torture = TRUE)
shine(x) ## interactive visualisation
```

![lineprof visualisation](https://camo.githubusercontent.com/13db9bc3ece496863d05c528c1d729d1f630247c/687474703a2f2f692e696d6775722e636f6d2f6e53437471734d2e706e67)

`lineprof` displays five variables for each line of code:

- `t`: the amount of time spent on that line (in seconds)
- `r`, `a`: the amount of memory released and allocated (in
  megabytes). The assignment of memory release to a line of is not
  deterministic because it occurs only when `gc` is triggered.
- `d`: the number of duplicates

### See also

The [`profr`](http://cran.r-project.org/web/packages/profr/index.html)
and
[`proftools`](http://cran.r-project.org/web/packages/proftools/index.html)
packages.

## Memory profiling

Memory usage using `tracemem` (requires to build `R` with `--enable-memory-profiling`)


```r
library("sequences")
(a <- new("DnaSeq", sequence = "GCATCAGCAGCT"))
```

```
## Object of class DnaSeq 
##  Id:  
##  Length: 12 
##  Alphabet: A C G T 
##  Sequence: GCATCAGCAGCT
```

```r
tracemem(a)
```

```
## [1] "<0x4601f98>"
```

```r
seq(a) <- "GATC"
```

```
## tracemem[0x4601f98 -> 0x425bc08]: eval eval withVisible withCallingHandlers doTryCatch tryCatchOne tryCatchList tryCatch try handle evaluate_call evaluate in_dir block_exec call_block process_group.block process_group withCallingHandlers process_file knit 
## tracemem[0x425bc08 -> 0x41f3ba0]: seq<- seq<- eval eval withVisible withCallingHandlers doTryCatch tryCatchOne tryCatchList tryCatch try handle evaluate_call evaluate in_dir block_exec call_block process_group.block process_group withCallingHandlers process_file knit
```

The illusion of copying


```r
x <- 1:10
tracemem(x)
```

```
## [1] "<0x3d29538>"
```

```r
y <- x 
## 2 'copies' of x trigger a real copy
x[1] <- 1L
```

```
## tracemem[0x3d29538 -> 0x3cb9e50]: eval eval withVisible withCallingHandlers doTryCatch tryCatchOne tryCatchList tryCatch try handle evaluate_call evaluate in_dir block_exec call_block process_group.block process_group withCallingHandlers process_file knit
```

```r
## Only one copy of x
x[1] <- 2L
```

## Object sizes

Approximate object's size


```r
object.size(rnorm(10000))
```

```
## 80040 bytes
```

```r
print(object.size(rnorm(10000)), units="auto")
```

```
## 78.2 Kb
```

```r
print(object.size(rnorm(1000000)), units="auto")
```

```
## 7.6 Mb
```

## Byte compiling `R` code

The `compile` package - the `cmpfun` function compiles the body of a
closure and returns a new closure with the same formals and the body
replaced by the compiled body expression.


```r
library("compiler")
readFastaCmp <- cmpfun(readFasta)
f <- dir(system.file("extdata",package="sequences"),
         pattern="fasta", full.names=TRUE)

library("microbenchmark")
microbenchmark(readFasta(f), readFastaCmp(f), times = 1e2)
```

```
## Unit: microseconds
##             expr     min       lq     mean  median      uq      max neval
##     readFasta(f) 603.590 618.1235 656.9715 637.179 673.671  946.715   100
##  readFastaCmp(f) 600.057 613.5380 703.0776 626.811 687.359 3110.835   100
##  cld
##    a
##    a
```
Fibonacci example


```r
fibR <- function(n) {
    if (n == 0) return(0)
    if (n == 1) return(1)
    return( fibR(n - 1) + fibR(n - 2))
}
fibR(10)
```

```
## [1] 55
```


```r
library("compiler")
fibRcmp <- cmpfun(fibR)
fibRcmp(10)
```

```
## [1] 55
```




```r
## a C++ implementation (see later)
fibC(10)
```

```
## [1] 55
```


```r
microbenchmark(fibR(10), fibRcmp(10), fibC(10), times = 1e2)
```

```
## Unit: microseconds
##         expr     min       lq      mean   median       uq      max neval
##     fibR(10) 152.605 175.9425 193.23615 178.7570 188.5425 1146.152   100
##  fibRcmp(10) 160.187 173.4690 182.94544 178.2695 185.3835  240.485   100
##     fibC(10)   1.950   2.6965   3.33669   3.1990   3.8300    8.476   100
##  cld
##    b
##    b
##   a
```



## Calling foreign languages

- R is getting too slow or is not doing well in terms of memory
  management.
- Implement the heavy stuff in `C`, `C++` (http://www.rcpp.org/),
  `Fortran` or `Java` (http://www.rforge.net/rJava/).


### Other scripting languages
- R/Perl: http://www.omegahat.org/RSPerl/ and 
- R/Python: http://www.omegahat.org/RSPython/ bidirectional interfaces

- There is also the `system()` function for direct access to OS functions

See C/C++ slides for details.

## Parallel execution

See [topic section](https://github.com/lgatto/rbc/tree/2014-11-06-Zurich/R-parallel).

## Environments

- An environment is a *isolated* collection of named objects, and a
  pointer to an enclosing environment.

- When calling a function, for example, its code is executed in an new
  environment; the function variables are local to that environment
  and distinct to those in your workspace (the `.GlobalEnv`).


```r
x <- 1
environment()
```

```
## <environment: R_GlobalEnv>
```

```r
f <- function() { print(environment()); x <- 2 }
f()
```

```
## <environment: 0x3f58cb8>
```

```r
x
```

```
## [1] 1
```


- We can create and populate new environments:

```r
e <- new.env()
e
```

```
## <environment: 0x3c98228>
```

```r
assign("x", value = 1, envir = e)
ls(envir = e)
```

```
## [1] "x"
```

```r
get("x", envir = e)
```

```
## [1] 1
```

- As well as lock/unlock bindings (with functions `lockBinding` and
  `unlockBonding`) or full environments (with `lockEnvironment`).

### Environments as function arguments

When passing an environment as function argument, it is **not**
copied: all its values are accessible within the function and can be
persistently modified.


```r
e <- new.env()
e$x <- 1
f <- function(myenv) myenv$x <- 2
f(e)
e$x
```

```
## [1] 2
```

This is used in the `eSet` et al. microarray data structures to store the expression data.

## Big data

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
- ...




