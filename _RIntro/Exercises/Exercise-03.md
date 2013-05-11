# How to store microarray data?

## The expression data

The most natural way to store the expression data is with a `matrix`,
where genes are represented along the rows and samples along the columns.
All values are of the same type (numerical values). 
Below we generate a matrix of random (normally distributed data),
representing 50 genes measured for 4 samples.


```r
expdata <- matrix(rnorm(200), nrow = 50, ncol = 4)
head(expdata)
```

```
##         [,1]     [,2]    [,3]     [,4]
## [1,]  0.3148 -0.79262  1.7571 -0.07668
## [2,]  0.4207 -1.58493 -0.8764 -0.07054
## [3,] -0.9228 -2.34246 -0.8919  0.03740
## [4,]  1.3592 -0.04233  0.6697  0.36964
## [5,]  1.3146 -1.62385  1.9986  0.42057
## [6,] -0.3190  0.75425  0.4131  0.64756
```


Note that this also support missing data encoded as `NA`.


```r
expdatana <- expdata  ## a copy of expdata
expdatana[2, 2] <- NA
head(expdatana)
```

```
##         [,1]     [,2]    [,3]     [,4]
## [1,]  0.3148 -0.79262  1.7571 -0.07668
## [2,]  0.4207       NA -0.8764 -0.07054
## [3,] -0.9228 -2.34246 -0.8919  0.03740
## [4,]  1.3592 -0.04233  0.6697  0.36964
## [5,]  1.3146 -1.62385  1.9986  0.42057
## [6,] -0.3190  0.75425  0.4131  0.64756
```

```r
class(expdatana)
```

```
## [1] "matrix"
```

```r
mode(expdatana)
```

```
## [1] "numeric"
```


We also want to set gene (row) and sample (column) names.


```r
dimnames(expdata) <-
  list(features = paste0("gene", 1:nrow(expdata)),
       samples = paste0("sample", 1:ncol(expdata)))
```


## The meta data

We provide addition information on the genes and samples using two data frames. 
We have to make sure that the respective dimensions of the expression data 
and gene/sample meta data match.


```r
smdata <- data.frame(feature = colnames(expdata),
                     group = c("ctrl", "ctrl",
                       "cond1", "cond1"),
                     replicate = rep(1:2, 2))
smdata
```

```
##   feature group replicate
## 1 sample1  ctrl         1
## 2 sample2  ctrl         2
## 3 sample3 cond1         1
## 4 sample4 cond1         2
```

```r
class(smdata)
```

```
## [1] "data.frame"
```

```r
nrow(smdata)
```

```
## [1] 4
```

```r
ncol(expdata)
```

```
## [1] 4
```

```r
nrow(smdata) == ncol(expdata)
```

```
## [1] TRUE
```



```r
fmdata <- data.frame(feature = rownames(expdata),                     
                     description = paste("Important gene", rownames(expdata)))
fmdata
```

```
##    feature           description
## 1    gene1  Important gene gene1
## 2    gene2  Important gene gene2
## 3    gene3  Important gene gene3
## 4    gene4  Important gene gene4
## 5    gene5  Important gene gene5
## 6    gene6  Important gene gene6
## 7    gene7  Important gene gene7
## 8    gene8  Important gene gene8
## 9    gene9  Important gene gene9
## 10  gene10 Important gene gene10
## 11  gene11 Important gene gene11
## 12  gene12 Important gene gene12
## 13  gene13 Important gene gene13
## 14  gene14 Important gene gene14
## 15  gene15 Important gene gene15
## 16  gene16 Important gene gene16
## 17  gene17 Important gene gene17
## 18  gene18 Important gene gene18
## 19  gene19 Important gene gene19
## 20  gene20 Important gene gene20
## 21  gene21 Important gene gene21
## 22  gene22 Important gene gene22
## 23  gene23 Important gene gene23
## 24  gene24 Important gene gene24
## 25  gene25 Important gene gene25
## 26  gene26 Important gene gene26
## 27  gene27 Important gene gene27
## 28  gene28 Important gene gene28
## 29  gene29 Important gene gene29
## 30  gene30 Important gene gene30
## 31  gene31 Important gene gene31
## 32  gene32 Important gene gene32
## 33  gene33 Important gene gene33
## 34  gene34 Important gene gene34
## 35  gene35 Important gene gene35
## 36  gene36 Important gene gene36
## 37  gene37 Important gene gene37
## 38  gene38 Important gene gene38
## 39  gene39 Important gene gene39
## 40  gene40 Important gene gene40
## 41  gene41 Important gene gene41
## 42  gene42 Important gene gene42
## 43  gene43 Important gene gene43
## 44  gene44 Important gene gene44
## 45  gene45 Important gene gene45
## 46  gene46 Important gene gene46
## 47  gene47 Important gene gene47
## 48  gene48 Important gene gene48
## 49  gene49 Important gene gene49
## 50  gene50 Important gene gene50
```

```r
nrow(fmdata)
```

```
## [1] 50
```

```r
nrow(expdata)
```

```
## [1] 50
```

```r
nrow(fmdata) == nrow(expdata)
```

```
## [1] TRUE
```


## A complete representation

We now combine the respective objects into a list to keep the 
expression data and its description together. 
We then print the structure of the list as a summary to ensure that 
we have the expected representation.



```r
marray <- list(
  expression = expdata,
  featuremeta = fmdata,
  samplemeta = smdata)
str(marray)
```

```
## List of 3
##  $ expression : num [1:50, 1:4] 0.315 0.421 -0.923 1.359 1.315 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ features: chr [1:50] "gene1" "gene2" "gene3" "gene4" ...
##   .. ..$ samples : chr [1:4] "sample1" "sample2" "sample3" "sample4"
##  $ featuremeta:'data.frame':	50 obs. of  2 variables:
##   ..$ feature    : Factor w/ 50 levels "gene1","gene10",..: 1 12 23 34 45 47 48 49 50 2 ...
##   ..$ description: Factor w/ 50 levels "Important gene gene1",..: 1 12 23 34 45 47 48 49 50 2 ...
##  $ samplemeta :'data.frame':	4 obs. of  3 variables:
##   ..$ feature  : Factor w/ 4 levels "sample1","sample2",..: 1 2 3 4
##   ..$ group    : Factor w/ 2 levels "cond1","ctrl": 2 2 1 1
##   ..$ replicate: int [1:4] 1 2 1 2
```

