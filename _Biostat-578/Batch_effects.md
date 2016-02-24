# Bioinformatics for Big Omics Data: The danger of batch effects!
Raphael Gottardo  
February 2, 2014  


## Setting up some options

Let's first turn on the cache for increased performance and improved styling

```r
# Set some global knitr options
library("knitr")
opts_chunk$set(tidy=TRUE, tidy.opts=list(blank=FALSE, width.cutoff=60), cache=TRUE, messages=FALSE)
```



## Reading

Before we start, you should read the following papers:

1. Leek, J. T. & Storey, J. D. Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis. PLoS Genet 3, e161 (2007).
2. Leek, J. T. & Storey, J. D. A general framework for multiple testing dependence. Proc. Natl. Acad. Sci. U.S.A. 105, 18718–18723 (2008).
3. Leek, J. T. et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nature Reviews Genetics 11, 733-739 (2010).
4. Johnson, W. E., Li, C. & Rabinovic, A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics 8, 118-127 (2007).
5. Gagnon-Bartsch, J. A. & Speed, T. P. Using control genes to correct for unwanted variation in microarray data. Biostatistics 13, 539-552 (2012).

## Motivation

Batch effects are technical sources of variation that have been added to the samples during handling.
Example of batch variables: lot number, technician, instrument settings, etc. 

If not adjusted for, these batch variables can have a substantial effects on downstream analysis


## Normalization and batch effects

Unfortunately, normalization will not always correct for batch effects. Technical variation due to batch effects might only affect a subset of the genes.

<img src="http://www.nature.com/nrg/journal/v11/n10/images/nrg2825-f1.jpg" width=300>

"For a published bladder cancer microarray data set obtained using an Affymetrix platform, we obtained the raw data for only the normal samples. Here, green and orange represent two different processing dates."

## Adjusting for batch effects

**Two scenarios:**

1. You have information about the batch variable
    - Use your batch effect as a covariate in your analysis (e.g. limma)
2. You suspect a batch effect, but you don't know where it is coming from
    - The batch effect needs to be estimated first and then corrected for, by adding the estimated variables as co-variates


## Singular value decomposition

Let $X$ be a matrix of size $m\times n$ ($m \ge n$) and rank $r\le n$
then we can decompose $X$ as 

$$X=USV^T$$

- U is the matrix of left singular vectors (eigenassays)
- V is the matrix of right singular vectors (eigengenes)
- S is the matrix of singular values (eigenvalues)

$U^TU=VV^T=I$ (orthogonal vectors)


$S=diag(s_1, \dots, s_n)$ where $s_l\ge 0$ and $s_{r+1}=\dots=s_n=0$

$X_i=\sum_j u_{ij}s_j\mathbf{v}_j$, which can be interpreted as a change of coordinate

## Relationship to principale component analysis

$X=USV^T$, and we have $X^TX=VSU^TUSV^T=VS^2V^T$

What happens if the rows of X are scaled?


## Surrogate variable analysis

Let $X_{m\times n}=(x_1,..,x_m)^T$ be the matrix of normalized expression values, with $n$ arrays and $m$ genes. 
Let $y=(y_1,..,y_n)^T$ be a vector of length $n$ representing the primary variable of interest (e.g covariates, vector of unknown coefficients). Without loss of generality model $x_{ij}=\mu_i+f_i( y_j) + e_{ij}$, where $\mu_i$ is the baseline level of expression, $f_i(y_j)=\mathbb{E}(x_{ij} | y_j)-\mu_i$ gives the relationship between measured variable of interest and gene $i$, and $e_{ij}$ is random noise with mean zero.

Suppose there are $L$ biologically meaningful unmodeled factors, such as age, environmental exposure, genotype,
etc. Let $g_l = (g_{l1},...,g_{ln})$ be an _arbitrarily complicated function_ of
the lth factor across all $n$ arrays, for $l=1,2,...,L$. Our model becomes:

$$x_{ij}=\mu_i + f_i(y_j) +\sum_{l=1}^L \gamma_{l_i}g_{l_j} + e^*_{ij}$$

and if factor $l$ does not influence the expression of gene $i$, we have $\gamma_{l_i}=0$.

## Surrogate variable analysis

In practice it is impossible to estimate $\sum_{l=1}^L \gamma_{l_i}g_{l_j}$, so Leek and Storey propose to use singular value decomposition to approximate the matrix $(\sum_{l=1}^L \gamma_{l_i}g_{l_j})_{ij}$ by its singular value decomposition. Computationally, this is done in two steps:


1. Detect unmodeled factors
2. Construct surrogate variables

## Detect unmodel factor
  
The main idea is as follows:

- Compute the residual matrix $r_{ij} = x_{ij}- \hat{\mu}_i - \hat{f}_i(y_j)$
- Perform the SVD of $R=(r_{ij})$
- Permute the rows of the matrix $R$ to obtain $R^*$. Regress  $r^*_{ij} = x_{ij}- \hat{\mu}_i - \hat{f}_i(y_j)$ to get residual matrix $R_0$, and perform the SVD of $R_0$. Repeat this many times to generate a null distribution for the residuals, given that $y$ is accounted for. 
- Compare the observed eigenvalues to those generated from the null distribution to obtain significance p-values
- Record the $K$ significant variables

## Construct surrogate variables 

1. Compute the residual matrix $r_{ij} = x_{ij} - \hat{\mu}_i - \hat{f}_i(y_j)$
2. Perform the SVD of $R=(r_{ij})$

Let $e_k=(e_{k1},...,e_{kn})^T$ be the $k$-th column of $V$ (for $k=1,...,n$). These $e_k$ are the residual eigengenes and represent orthogonal residual signals independent of the signal due to the primary variable.

2. Regress $e_k$ on the $x_i$ to access the significance of the $k$-th factor on gene $i$
3. Use the selected genes to form a reduced expression matrix and repeat 1. The estimated factor will form the basis for the surrogate variables
4. In any subsequent analysis include these factors in your model

## Using the sva package

We're going to look at the dataset used in:

Nakaya, H. I., Wrammert, J., Lee, E. K., Racioppi, L., Marie-Kunze, S., Haining, W. N., et al. (2011). Systems biology of vaccination for seasonal influenza in humans. Nature Immunology, 12(8), 786795. doi:10.1038/ni.2067


```r
library(GEOquery)
```

```
## Loading required package: Biobase
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Setting options('download.file.method.GEOquery'='auto')
```

```r
# Download the mapping information and processed data
gds <- getGEO("GSE29619", destdir = "Data/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29619/matrix/
## Found 3 file(s)
## GSE29619-GPL13158_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL13158_series_matrix.txt.gz
## Using locally cached version of GPL13158 found here:
## Data/GEO//GPL13158.soft 
## GSE29619-GPL3921_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL3921_series_matrix.txt.gz
## Using locally cached version of GPL3921 found here:
## Data/GEO//GPL3921.soft 
## GSE29619-GPL570_series_matrix.txt.gz
## Using locally cached version: Data/GEO//GSE29619-GPL570_series_matrix.txt.gz
## Using locally cached version of GPL570 found here:
## Data/GEO//GPL570.soft
```

```r
# main serie #gds[[1]] = LAIV/TIV 0809, gds[[2]] = FACS,
# gds[[3]] = TIV 0708
```

but before we can use this, we need to clean up the pData a bit (see code in .Rpres file).



## Using the sva package 

Let us estimate the (surrogate) factors as follows:

```r
library(sva)
```

```
## Loading required package: corpcor
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.8-4. For overview type 'help("mgcv-package")'.
```

```r
TIV_08 <- gds_new[[1]][, grepl("2008-TIV", pData(gds_new[[1]])$description)]
mm_TIV_08 <- model.matrix(~ptid + time, TIV_08)
mm0_TIV_08 <- model.matrix(~ptid, TIV_08)
# Estimate the surrogate variables
sv_TIV_08 <- sva(exprs(TIV_08), mm_TIV_08, mod0 = mm0_TIV_08)
```

```
## Number of significant surrogate variables is:  12 
## Iteration (out of 5 ):1  2  3  4  5
```

## Using the limma with estimated SVs

Then we can use these variables in limma, as follows:

```r
library(limma)
```

```
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
# Add the surrogate variables to the design matrix
modSv <- cbind(mm_TIV_08, sv_TIV_08$sv)
# Use the new design matrix
fit_TIV_08 <- lmFit(TIV_08, modSv)
ebay_TIV_08 <- eBayes(fit_TIV_08)
topT7_sv <- topTable(ebay_TIV_08, coef = "timeD7", number = Inf)
# Compare to the old analysis
fit_TIV_08 <- lmFit(TIV_08, mm_TIV_08)
ebay_TIV_08 <- eBayes(fit_TIV_08)
topT7 <- topTable(ebay_TIV_08, coef = "timeD7", number = Inf)
```

## Using the limma with estimated SVs

The result from the adjusted analysis:


```r
topT7_sv[1:10, c("ID", "Gene.Symbol")]
```

```
##                            ID
## 211430_PM_s_at 211430_PM_s_at
## 232991_PM_at     232991_PM_at
## 219276_PM_x_at 219276_PM_x_at
## 239401_PM_at     239401_PM_at
## 239637_PM_at     239637_PM_at
## 227999_PM_at     227999_PM_at
## 226481_PM_at     226481_PM_at
## 216801_PM_at     216801_PM_at
## 241402_PM_at     241402_PM_at
## 235180_PM_at     235180_PM_at
##                                                          Gene.Symbol
## 211430_PM_s_at IGH@ /// IGHG1 /// IGHM /// IGHV4-31 /// LOC100290146
## 232991_PM_at                                                        
## 219276_PM_x_at                                               C9orf82
## 239401_PM_at                                                        
## 239637_PM_at                                                        
## 227999_PM_at                                                  PWWP2B
## 226481_PM_at                                                   VPRBP
## 216801_PM_at                                                        
## 241402_PM_at                                                  TSEN54
## 235180_PM_at                                                    STYX
```

## Using the limma with estimated SVs

The result from the un-adjusted analysis:


```r
topT7[1:10, c("ID", "Gene.Symbol")]
```

```
##                              ID
## 211430_PM_s_at   211430_PM_s_at
## 241824_PM_at       241824_PM_at
## 216576_PM_x_at   216576_PM_x_at
## 1559018_PM_at     1559018_PM_at
## 216207_PM_x_at   216207_PM_x_at
## 214669_PM_x_at   214669_PM_x_at
## 1568768_PM_s_at 1568768_PM_s_at
## 215946_PM_x_at   215946_PM_x_at
## 211645_PM_x_at   211645_PM_x_at
## 217157_PM_x_at   217157_PM_x_at
##                                                           Gene.Symbol
## 211430_PM_s_at  IGH@ /// IGHG1 /// IGHM /// IGHV4-31 /// LOC100290146
## 241824_PM_at                                                         
## 216576_PM_x_at              IGK@ /// IGKC /// LOC652493 /// LOC652694
## 1559018_PM_at                                                   PTPRE
## 216207_PM_x_at                                              IGKV1D-13
## 214669_PM_x_at            IGK@ /// IGKC /// IGKV3-20 /// LOC100291682
## 1568768_PM_s_at                                          LOC100302650
## 215946_PM_x_at                           IGLL1 /// IGLL3 /// LOC91316
## 211645_PM_x_at                                                       
## 217157_PM_x_at                            IGK@ /// IGKC /// LOC652493
```

What do you think?
