# R and Bioconductor for proteomics

## Bioconductor

[Bioconductor](http://www.bioconductor.org/) provides tools for the
analysis and comprehension of high-throughput biology (genomics,
proteomics, metabolomics, flow cytometry, ...)  data. Bioconductor
uses the R statistical programming language, and is open source and
open development. It has two releases each year, 1295 software
packages, and an active user community.

In addition to *software* package, it also has *experiment* and
*annotation* packages.

## Installing Bioconductor packages

First time


```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RforProteomics")
```

Then


```r
library("BiocInstaller")
biocLite("RforProteomics")
```

To install all Bioconductor software that will be used throughout this
workshop:


```r
biocLite(c("RforProteomics", "pRolocGUI"), dependencies = TRUE)
```

See the detailed session information in the [wrap up](../wrapup.md)
section for packages and versions used to generate this material.

## Support

- the [Bioconductor support forum](http://support.bioconductor.org/)

Some packages have GitHub pages and use GitHub issues. That would be a
good place to report bugs. But the preferred place to ask questions is
the support forum.

## Proteomics and mass spectrometry

The easiest way to browser and search Bioconductor packages is to
navigate [task
views](http://www.bioconductor.org/packages/release/BiocViews.html#___Software).

Particularly relevant for this course are the [mass
spectrometry](http://www.bioconductor.org/packages/release/BiocViews.html#___MassSpectrometry),
[mass spectrometry
data](http://www.bioconductor.org/packages/release/BiocViews.html#___MassSpectrometryData)
and
[proteomics](http://www.bioconductor.org/packages/release/BiocViews.html#___Proteomics)
task views.

## The `RforProteomics` package

The
[`RforProteomics`](http://www.bioconductor.org/packages/release/data/experiment/html/RforProteomics.html)
package contains code to illustrate the *Using R and Bioconductor for
proteomics data analysis* and *Visualisation of proteomics data using
R and Bioconductor* papers. Two vignettes describe the code and data
needed to reproduce all the examples and figures.





```r
library("RforProteomics")
```

```
## 
## This is the 'RforProteomics' version 1.11.2.
## 
##   To get started, visit
##     http://lgatto.github.com/RforProteomics/
## 
##   or, in R, open package vignettes by typing
##     RforProteomics() # R/Bioc for proteomics overview
##     RProtVis()       # R/Bioc for proteomics visualisation
## 
##   For a full list of available documents:
##     vignette(package='RforProteomics')
```

Package vignettes are overview package documentations. They generally
provide a complete use case demonstrating the package's
functionality. All Bioconductor software packages have vignettes in
addition to all function manuals.

## Package landing pages

Each Bioconductor package has an official page on the Bioconductor website:

* `http://www.bioconductor.org/packages/packageName`

For example

* http://www.bioconductor.org/packages/RforProteomics
* http://www.bioconductor.org/packages/MSnbase

These pages summarise general information about the package and
provides links to its vignettes.

## Proteomics/MS data structures


|Data type      |File format          |Data structure               |Package           |
|:--------------|:--------------------|:----------------------------|:-----------------|
|Raw            |mzXML or mzML        |mzRpwiz or mzRramp           |mzR               |
|Raw            |mzXML or mzML        |list of MassSpectrum objects |MALDIquantForeign |
|Raw            |mzXML or mzML        |MSnExp                       |MSnbase           |
|Identification |mzIdentML            |mzRident                     |mzR               |
|Identification |mzIdentML            |mzID                         |mzID              |
|Quantitative   |mzTab                |MSnSet                       |MSnbase           |
|Peak lists     |mgf                  |MSnExp                       |MSnbase           |
|Imaging        |imzML or Analyze 7.5 |MSImageSet                   |Cardinal          |
|Imaging        |imzML or Analyze 7.5 |list of MassSpectrum objects |MALDIquantForeign |

## References

Gatto L. and Christoforou A. *Using R and Bioconductor for proteomics
data analysis*, Biochim Biophys Acta - Proteins and
Proteomics, 2013. [PMID:23692960](https://www.ncbi.nlm.nih.gov/pubmed/23692960)([preprint](https://arxiv.org/abs/1305.6559))

Gatto L, Breckels LM, Naake T, Gibb S. *Visualisation of proteomics
data using R and Bioconductor*. Proteomics. 2015 Feb 18. doi:
10.1002/pmic.201400392. [PMID:25690415](http://www.ncbi.nlm.nih.gov/pubmed/25690415).



