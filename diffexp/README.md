# Identifying differentially expressed proteins




## Running a `t-test` in R

Usig the `t.test` function:

```
t.test(x, y = NULL,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95, ...)
```

We will focus on two sample unpaires t-test, assuming unequal
variances, as this is the most common scenario in proteomics. Using a
**paired test** when appropriate is essential, as it will
substantially increase your test power.

We are going to use the `rnorm` function in this an the next section
to quickly genreate normally distributed data. Its inputs are 

- `n`: the number of data points to be generated;
- `mean`: the mean of the normal distribution to draw the data from
  (default is 0);
- `sd`: the standard deviation of the normal distribution to draw the
  data from (default is 1).

### Exercise 

* Generate 200 numbers drawn from a normal distribution of mean 0 and
  standard deviation 1. Verify that the parameters of the randomly
  data are correct. What figure would you use to visualise such data?
  
* Same as above for a normal distribution of mean 2 and standard
  deviation 0.5.

* Do you get the same values as your neighbour?

Let's now apply a t-test on two sets of values drawn from identical
and different distributions:


```r
t1 <- t.test(rnorm(5), rnorm(5))
t1
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  rnorm(5) and rnorm(5)
## t = 0.5154, df = 5.9441, p-value = 0.6249
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -1.292701  1.980590
## sample estimates:
##   mean of x   mean of y 
##  0.04179283 -0.30215135
```


```r
t2 <- t.test(rnorm(5), rnorm(5, mean = 4))
t2
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  rnorm(5) and rnorm(5, mean = 4)
## t = -11.458, df = 7.5818, p-value = 4.678e-06
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -5.288970 -3.502512
## sample estimates:
##  mean of x  mean of y 
## -0.2739408  4.1217999
```

What we see above is a pretty output that is convenient to visualise
interactively. The output of the `t.test` is an object of class `r
class(t2)`, which contains statistic, parameter, p.value, conf.int, estimate, null.value, alternative, method, data.name.



## Multiple testing



## Visualising results

Volcano plot


## Moderated t-tests: `limma`

## Count data

## Other packages

* *[MSstats](http://bioconductor.org/packages/MSstats)* for various statistical analyses
* Isobaric tagging (iTRAQ and TMT): *[isobar](http://bioconductor.org/packages/isobar)*
* Label-free: *[aLFQ](http://cran.fhcrc.org/web/packages/aLFQ/index.html)* and *[protiq](http://cran.fhcrc.org/web/packages/protiq/index.html)*

