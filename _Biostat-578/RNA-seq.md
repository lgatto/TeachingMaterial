# Bioinformatics for Big Omics Data: RNA-seq data analysis
Raphael Gottardo  
February 11, 2015  


## Setting up some options

Let's first turn on the cache for increased performance and improved styling

```r
# Set some global knitr options
library("knitr")
opts_chunk$set(tidy=TRUE, tidy.opts=list(blank=FALSE, width.cutoff=60), cache=TRUE, messages=FALSE)
```


We will be using the following packages

```r
library(data.table)
library(ggplot2)
library(limma)
library(edgeR)
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
## The following object is masked from 'package:limma':
## 
##     plotMA
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     Filter, Find, Map, Position, Reduce, anyDuplicated, append,
##     as.data.frame, as.vector, cbind, colnames, do.call,
##     duplicated, eval, evalq, get, intersect, is.unsorted, lapply,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, rank, rbind, rep.int, rownames, sapply, setdiff,
##     sort, table, tapply, union, unique, unlist, unsplit
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Setting options('download.file.method.GEOquery'='curl')
```


## Outline

Here we will discuss RNA-seq data analysis including normalization and differential expression. You should read the following papers:

1. Marioni, J. C., Mason, C. E., Mane, S. M., Stephens, M. & Gilad, Y. RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays. Genome Res. 18, 1509-1517 (2008).

2. Robinson, M. D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol. 11, R25 (2010).

3. Anders, S. et al. Count-based differential expression analysis of RNA sequencing data using R and Bioconductor. Nat Protoc 8, 1765-1786 (2013).

4. Lund, S. P., Nettleton, D., McCarthy, D. J. & Smyth, G. K. Detecting differential expression in RNA-sequence data using quasilikelihood with shrunken dispersion estimates. Stat Appl Genet Mol Biol 11, (2012).

5. Anders, S. & Huber, W. Differential expression analysis for sequence count data. Genome Biol. 11, R106 (2010).

7. Law, C. W., Chen, Y., Shi, W. & Smyth, G. K. Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 15, R29 (2014).

## RNA-seq workflow


<img src="Images/RNA-seq-workflow.png" width=500>

## Modeling RNA-seq data


As opposed to microarrays, sequencing leads to count data. In our case, we will try to model counts over genes (or possibly genomic intervals). Let's denote by $Y_{ij}$ the number of reads that fall withing gene $i$ in sample $j$. Then a possible discrete model is the Poisson distribution

$$
Y_{ijk}=\mathrm{Poisson}(w_{ij}N_{jk})
$$

where $w_{ij}$ is the normalized rate of emission of gene $i$ in sample $j$, and $N_{jk}$ is the total read count in lane $k$ for sample $j$. $N_{jk}$ basically accounts for sequencing depth variability from one lane to the next. This is basically the model proposed by Marioni et al. (2010). Note that we have the following constraint $\sum_i w_{ij}=1$.


## Modeling RNA-seq data 


The disavantage of the Poisson distribution is that it depends on one parameter only, which can lead to a lack of fit (over-dispersion). Many groups have proposed the negative binomial as an alternative model:

In particular, Robinson and Smyth (2007,2008) proposed the following model:

$$
Y_{ijk}=\mathrm{NB}(w_{ij}N_{jk}, \phi)
$$

where $w_{ij}$, $N_{jk}$ are as defined previously and $\phi$ is the dispersion parameter. Robinson and Smyth use the following parametrization: 

$$\mathbb{E}(Y_{ijk})=\mu_{ijk}\quad \mathrm{and}\quad \mathrm{Var}(Y_{ijk})=\mu_{ijk}(1+\mu_{ijk}\phi)$$

where $\mu_{ijk}=w_{ij}N_{jk}$. In this parametrization the poisson distribution can be seen as a special case when $\phi=0$.

The genewise dispersion parameters are estimated by conditional maximum likelihood, conditioning on the total count for that gene (Smyth and Verbyla, 1996). An empirical Bayes procedure is used to shrink the dispersions towards a consensus value, effectively borrowing information between genes (Robinson and Smyth, 2007).

## Estimating the dispersion parameter


As in `limma`, we would like to use an EB approach to shrink the gene-wise dispersion parameters towards a common value. Unfortunately, the NB distribution does not have a conjugate prior. Robinson and Smyth (2007) proposed to use weighted likelihood. Instead of maximizing the likelihood, they propose to maximize the following weighted likelihood

$$
WL(\phi_g)=l_g(\phi_g)+\alpha l_C(\phi_g)
$$

where $\alpha$ is the weight given to the common likelihood (i.e. $\phi$ constant across genes). Then the problem becomes the choice of the appropriate $\alpha$ value. R&S propose some strategies using some EB arguments. 

## edgeR

The model described above is implemented in `edgeR`, which also provides a generalized linear model interface for modeling the mean expression value as a function of explanatory variables. This is very similar to the `limma` framework. Then inference is done using a likelihood ratio-test comparing the alternative model to the null model. 


## DESeq

`DEseq` is another popular package for RNA-seq analysis, which also utilizes an NB model. 
In `edgeR`, the mean variance relationship is defined as $\sigma^2=\mu+\alpha \mu^2$. DESeq differs from `edgeR` in its mean-variance relationship, and the way the dispersion parameters are estimated. As explained in Anders et al. (2013):

> edgeR moderates feature-level dispersion estimates toward a trended mean according to the dispersion-mean relationship. In contrast, DESeq takes the maximum of the individual dispersion estimates
and the dispersion-mean trend. In practice, this means DESeq is less powerful, whereas edgeR is more sensitive to outliers. Recent comparison studies have highlighted that no single method dominates another across all settings.

## Using a normal approximation?


An alternative to the Poisson and NB models would be to find a data transformation that would make the count data approximately normal. Law et al. (2014) propose to use the $\log_2(\mathrm{cpm}(\cdot))$ transformation. The cpm transformation accounts for sequencing depth variability while the $\log_2$ transformation makes the data more normal. However, as we've seen earlier, the mean-variance relationship is quadratic, and a log transformation is not going to remove this dependence. As a consequence, it would be innapropriate to use a normal linear model with constant variance (even gene-wise). 

## Mean-variance trend estimation via voom

Law et al. (2014) propose to estimate the mean-variance trend, and then to use the inverse of the estimated standard deviation for each observation as weight in LIMMA. 
This is done by the `voom` function in LIMMA. 

Basically, `voom` fits the linear model without weights and uses the residuals of the model to estimate the weights, which are then pass onto a weighted LIMMA call for linear modeling. 

Law et al. (2014) show that this approach: 

1. Control the type I error rate
2. Is powerful among the methods that control the type I error rate
3. Has good FDR control
4. Is faster!

The `voom`+`limma` approach can also be used for gene set analysis, which is difficult to do with count based methods. 

## RNA-seq example

Get the data for GSE45735 from GEO.


```r
# You should make sure the directory Data/GEO exists
datadir <- "Data/GEO/"
dir.create(file.path(datadir), showWarnings = FALSE, recursive = TRUE)
# Get the data from GEO
gd <- getGEO("GSE45735", destdir = datadir)
pd <- pData(gd[[1]])
getGEOSuppFiles("GSE45735", makeDirectory = FALSE, baseDir = datadir)
```

## RNA-seq example

Clean up the data.


```r
# Note the regular expression to grep file names
files <- list.files(path = datadir, pattern = "GSE45735_T.*.gz", 
    full.names = TRUE)
# Read in gzip-compressed, tab-delimited files
file_list <- lapply(files, read.table, sep = "\t", header = TRUE)
# Subset to only those rows where Gene contains only
# non-space characters This addresses problems with T14 file
# containing 28 invalid rows at end of file
file_list <- lapply(file_list, function(file_list) subset(file_list, 
    grepl("^[^[:space:]]+$", Gene)))
# Remove duplicated rows
file_list_unique <- lapply(file_list, function(x) {
    x <- x[!duplicated(x$Gene), ]
    x <- x[order(x$Gene), ]
    rownames(x) <- x$Gene
    x[, -1]
})
```

## RNA-seq example

Create a matrix from the intersection of all genes and clean up the pData.


```r
# Take the intersection of all genes
gene_list <- Reduce(intersect, lapply(file_list_unique, rownames))
file_list_unique <- lapply(file_list_unique, "[", gene_list, 
    )
matrix <- as.matrix(do.call(cbind, file_list_unique))
# Clean up the pData
pd_small <- pd[!grepl("T13_Day8", pd$title), ]
pd_small$Day <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title), 
    "_"), "[", 2)
pd_small$subject <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title), 
    "_"), "[", 1)
colnames(matrix) <- rownames(pd_small)
```


## RNA-seq example 

Note that raw data files for sequencing experiments are available from the SRA database, which can be queried using the SRAdb package:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("SRAdb")
```

The resulting files are usually very large!


## Using Voom


Let's first create an eSet we can use:

```r
# Note that I add one to the count
new_set <- ExpressionSet(assayData = matrix + 1)
pData(new_set) <- pd_small
```

we now need to set-up our design matrix to estimate our weights:


```r
design <- model.matrix(~subject + Day, new_set)
new_set_voom <- voom(new_set, design = design)
```


```r
lm <- lmFit(new_set_voom, design)
eb <- eBayes(lm)
# Look at the other time-points
topTable(eb, coef = "DayDay1", number = 5)
```

```
##            logFC  AveExpr        t      P.Value   adj.P.Val        B
## TRIM25 0.3183383 7.417747 6.613100 6.878706e-08 0.001509050 7.951518
## IRF1   0.6689637 7.822953 6.084187 3.764288e-07 0.004030395 6.244832
## IRF9   0.4553882 6.700297 5.845475 8.117054e-07 0.004030395 5.573664
## ASPHD2 0.5128938 4.020599 5.740800 1.136747e-06 0.004060609 5.426160
## STAT1  0.8919256 8.825765 5.840196 8.256134e-07 0.004030395 5.393302
```

## Using edgeR


Let's see how we would get setup for edgeR


```r
# We don't have the count matrix, so let's try to 'estimate'
# the counts
counts <- ceiling(matrix * 10)
dge_list <- DGEList(counts = counts, group = pd_small$Day)
dge_list <- calcNormFactors(dge_list)
design <- model.matrix(~Day + subject, pd_small)
dge_list <- estimateGLMCommonDisp(dge_list, design)
dge_list <- estimateGLMTrendedDisp(dge_list, design)
```

```
## Loading required package: splines
```

```r
dge_list <- estimateGLMTagwiseDisp(dge_list, design)
fit <- glmFit(dge_list, design)
lrt <- glmLRT(fit, coef = "DayDay1")
detags <- rownames(topTags(lrt, n = 20))
```

## Using edgeR


```r
plotSmear(lrt, de.tags = detags)
```

![](RNA-seq_files/figure-html/unnamed-chunk-11-1.png) 


## Normalization


As we have seen above, `edgeR` provides a function for estimating the normalizing factor. 
The purpose of this normalization strategy is to correct for sequencing depth variability across libraries/lanes. A natural way to achieve this would be simply scale the counts by the inverse of the sum of the counts across genes (i.e. the total number of reads/lane). However, genes with extreme expression values might bias this total read number estimate. Robinson & Oshlack propose a different approach for estimating the normalizing factor, Trimmed Mean of M-values. By default, edgeR uses the TMM normalization strategy.


## Normalization - TMM


<img src="Images/TMM-motivation.png" width=800>

## Normalization - TMM 


<img src="Images/TMM-formula.png" width=800>

A trimmed mean is the average after removing the upper and lower x% of the data. The TMM procedure is doubly trimmed, by trimming both the M and A values. By default, we trim the Mg values by 30% and the Ag values by 5%, but these settings can be tailored to a given experiment. 

## DEseq and DEXseq


DEseq is also another popular approach for differential expression, and I will leave it up to you to try this package. One major technical advantage of RNA-seq vs microarrays is that you actually get expression level (count) data at the exon level. As such, using RNA-seq it is possible to detect differential usage of exons. This is what the DEXseq package tries to do. DEXseq models binned (exon) counts instead of gene counts.

<img src="http://genome.cshlp.org/content/22/10/2008/F1.large.jpg" width=300>


## Summary


We've explored some of the computational challenges involved in the analysis of RNA-seq data. There exists many different methods, tools and software packages for RNA-seq. I encourage you to explore some of these tools/methods. Some of these have been proposed for your final project.
