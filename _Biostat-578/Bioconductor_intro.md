# Bioinformatics for Big Omics Data: Introduction to Bioconductor
Raphael Gottardo  
January 21, 2014  

## Setting up some options

Let's first turn on the cache for increased performance and improved styling

```r
# Set some global knitr options
library("knitr")
opts_chunk$set(tidy=TRUE, tidy.opts=list(blank=FALSE, width.cutoff=60), cache=TRUE, messages=FALSE)
```


## The Bioconductor project

- [Bioconductor](http://www.bioconductor.org) is an open source, open development software project to provide tools for the analysis and comprehension of high-throughput genomic data. It is based primarily on the R programming language.

- Most Bioconductor components are distributed as R packages. The functional scope of Bioconductor packages includes the analysis of DNA microarray, sequence, flow, SNP, and other data.

## Project Goals

The broad goals of the Bioconductor project are:

- To provide widespread access to a broad range of powerful statistical and graphical methods for the analysis of genomic data.
- To facilitate the inclusion of biological metadata in the analysis of genomic data, e.g. literature data from PubMed, annotation data from Entrez genes.
- To provide a common software platform that enables the rapid development and deployment of extensible, scalable, and interoperable software.
- To further scientific understanding by producing high-quality documentation and reproducible research.
- To train researchers on computational and statistical methods for the analysis of genomic data.

## Getting started

Before running this presentation, you will need to:

 - set the current directory in the console (go to session -> Set Working Directory -> To Source File Location)

- run the following commands which take up time
    - source("http://bioconductor.org/biocLite.R")
    - biocLite()
    - getSQLiteFile() Note: this will take some time
    - getGEOSuppFiles("GSE29617", makeDirectory = TRUE, baseDir = "./")
    - untar("./GSE29617/GSE29617_RAW.tar", exdir="./GSE29617", tar = Sys.getenv("TAR"))

Note that I have updated the code so that you may not have to do what's described above.


## Getting started


```r
source("http://bioconductor.org/biocLite.R")
```

```
## Bioconductor version 2.14 (BiocInstaller 1.14.3), ?biocLite for
##   help
## A newer version of Bioconductor is available for this version of
##   R, ?BiocUpgrade for help
```

```r
# Install all core packages and update all installed packages
biocLite()
```

```
## BioC_mirror: http://bioconductor.org
## Using Bioconductor version 2.14 (BiocInstaller 1.14.3), R version
##   3.1.2.
## Old packages: 'dplyr', 'robustbase'
```

You can also install specific packages


```r
biocLite(c("GEOmetadb", "GEOquery"))
```


## The Gene Expression Omnibus (GEO)

The [Gene Expression Omnibus](http://www.ncbi.nlm.nih.gov/geo/) is an international public repository that archives and freely distributes microarray, next-generation sequencing, and other forms of high-throughput functional genomics data submitted by the research community.

The three main goals of GEO are to:

- Provide a robust, versatile database in which to efficiently store high-throughput functional genomic data
- Offer simple submission procedures and formats that support complete and well-annotated data deposits from the research community
- Provide user-friendly mechanisms that allow users to query, locate, review and download studies and gene expression profiles of interest

## Getting data from GEO

Before getting data from GEO, we need to see what data we want. For that we can use the `GEOmetadb` package. 


```r
library(GEOmetadb)
```

Remember that packages in Bioconductor are well documented with a vignette that can be access as follows:


```r
vignette("GEOmetadb")
```
or if the package contains multiple vignettes or a vignette with a non-standard name

```r
browseVignettes(package = "GEOmetadb")
```

## Finding the right data in GEO

Zhu, Y., Davis, S., Stephens, R., Meltzer, P. S., & Chen, Y. (2008). GEOmetadb: powerful alternative search engine for the Gene Expression Omnibus. Bioinformatics (Oxford, England), 24(23), 2798–2800. doi:10.1093/bioinformatics/btn520

GEOmetadb uses a SQLite database to store all metadata associate with GEO.

```r
## This will download the entire database, so can be slow
if (!file.exists("GEOmetadb.sqlite")) {
    # Download database only if it's not done already
    getSQLiteFile()
}
```

## Finding the right data in GEO


```r
geo_con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListTables(geo_con)
```

```
##  [1] "gds"               "gds_subset"        "geoConvert"       
##  [4] "geodb_column_desc" "gpl"               "gse"              
##  [7] "gse_gpl"           "gse_gsm"           "gsm"              
## [10] "metaInfo"          "sMatrix"
```


```r
dbListFields(geo_con, "gse")
```

```
##  [1] "ID"                   "title"                "gse"                 
##  [4] "status"               "submission_date"      "last_update_date"    
##  [7] "pubmed_id"            "summary"              "type"                
## [10] "contributor"          "web_link"             "overall_design"      
## [13] "repeats"              "repeats_sample_list"  "variable"            
## [16] "variable_description" "contact"              "supplementary_file"
```

## Finding a study

The basic record types in GEO include Platforms (GPL), Samples (GSM), Series (GSE) and DataSets (GDS)


```r
dbGetQuery(geo_con, "SELECT gse.ID, gse.title, gse.gse FROM gse WHERE gse.pubmed_id='21743478';")
```

```
##      ID
## 1 26409
## 2 26410
## 3 26412
## 4 26413
## 5 26414
##                                                                                                          title
## 1                  Time Course of Young Adults Vaccinated with Influenza TIV Vaccine during 2007/08 Flu Season
## 2                 Time Course of Young Adults Vaccinated with Influenza LAIV Vaccine during 2008/09 Flu Season
## 3                  Time Course of Young Adults Vaccinated with Influenza TIV Vaccine during 2008/09 Flu Season
## 4 FACS-sorted cells from Young Adults Vaccinated with Influenza TIV or LAIV Vaccines during 2008/09 Flu Season
## 5                                              Systems biology of vaccination for seasonal influenza in humans
##        gse
## 1 GSE29614
## 2 GSE29615
## 3 GSE29617
## 4 GSE29618
## 5 GSE29619
```

## Finding a study

What samples were used?


```r
dbGetQuery(geo_con, "SELECT gse.gse, gsm.gsm, gsm.title FROM (gse JOIN gse_gsm ON gse.gse=gse_gsm.gse) j JOIN gsm ON j.gsm=gsm.gsm WHERE gse.pubmed_id='21743478' LIMIT 5;")
```

```
##    gse.gse   gsm.gsm                                     gsm.title
## 1 GSE29614 GSM733816 2007 TIV subject ID 12 at D0 post-vaccination
## 2 GSE29614 GSM733817 2007 TIV subject ID 12 at D3 post-vaccination
## 3 GSE29614 GSM733818 2007 TIV subject ID 12 at D7 post-vaccination
## 4 GSE29614 GSM733819 2007 TIV subject ID 16 at D0 post-vaccination
## 5 GSE29614 GSM733820 2007 TIV subject ID 16 at D3 post-vaccination
```
gse_gsm contains the gse number that is associated with the gsm number. 
j is the name of a table that is created by joining gse and ges_gsm. Then j is joined with table gsm. 

## Finding a study

What about raw data?


```r
res <- dbGetQuery(geo_con, "SELECT gsm.gsm, gsm.supplementary_file FROM (gse JOIN gse_gsm ON gse.gse=gse_gsm.gse) j JOIN gsm ON j.gsm=gsm.gsm WHERE gse.pubmed_id='21743478' LIMIT 5;")
head(res)
```

```
##     gsm.gsm
## 1 GSM733816
## 2 GSM733817
## 3 GSM733818
## 4 GSM733819
## 5 GSM733820
##                                                                                                                                                  gsm.supplementary_file
## 1 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733816/suppl/GSM733816.CEL.gz;\tftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733816/suppl/GSM733816.chp.gz
## 2 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733817/suppl/GSM733817.CEL.gz;\tftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733817/suppl/GSM733817.chp.gz
## 3 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733818/suppl/GSM733818.CEL.gz;\tftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733818/suppl/GSM733818.chp.gz
## 4 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733819/suppl/GSM733819.CEL.gz;\tftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733819/suppl/GSM733819.chp.gz
## 5 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733820/suppl/GSM733820.CEL.gz;\tftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733820/suppl/GSM733820.chp.gz
```
raw data is contained in the supplementary files, which are listed in the gsm file. 

## Finding specific data

To get list of manufacturers:

```r
library(data.table)
manu <- data.table(dbGetQuery(geo_con, "SELECT manufacturer FROM gpl"))
manu[, list(length(manufacturer)), by = manufacturer]
```

```
##                                                           manufacturer V1
##    1:                                                               NA  1
##    2:                                                       Affymetrix  1
##    3:                                          BD Biosciences Clontech  1
##    4:                                                    Sigma Genosys  1
##    5: Advanced Technology Center Microarray Facility (NIH/NCI/CCR/ATC)  1
##   ---                                                                    
## 2074:                                                Biovue Technology  1
## 2075:                                     NanoString, Seattle, WA, USA  1
## 2076:                                                         Homemade  1
## 2077:                                           Qiagen (SABiosciences)  1
## 2078:                                              Affymetrix Panomics  1
```

## Finding specific data

To get supplementary file names ending with cel.gz from only manufacturer Affymetrix

```r
res <- dbGetQuery(geo_con, "SELECT gpl.bioc_package, gsm.title, gsm.series_id, gsm.gpl, gsm.supplementary_file FROM gsm JOIN gpl ON gsm.gpl=gpl.gpl WHERE gpl.manufacturer='Affymetrix' AND gsm.supplementary_file like '%CEL.gz';")
head(res)
```

```
##   bioc_package               title series_id   gpl
## 1       hu6800          BM_CD34-1a    GSE500 GPL80
## 2       hu6800          BM_CD34-1b    GSE500 GPL80
## 3       hu6800           BM_CD34-2    GSE500 GPL80
## 4       hu6800        GPBMC_CD34-1    GSE500 GPL80
## 5       hu6800        GPBMC_CD34-2    GSE500 GPL80
## 6       rgu34a CNS-SC_Inj24h-3A-s2    GSE464 GPL85
##                                                            supplementary_file
## 1    ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM575/suppl/GSM575.cel.gz
## 2    ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM576/suppl/GSM576.cel.gz
## 3    ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM577/suppl/GSM577.cel.gz
## 4    ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM578/suppl/GSM578.cel.gz
## 5    ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM579/suppl/GSM579.cel.gz
## 6 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1nnn/GSM1136/suppl/GSM1136.CEL.gz
```

## Finding specific data

From previous table:

- bioc_package = bioconductor package
- hu6800 = Affymetrix HuGeneFL Genome Array annotation data (chip hu6800) 
- rgu34a = Affymetrix Rat Genome U34 Set annotation data (chip rgu34a)
- title = data set title or study title

For example BM_CD34-1a = bone marrow flow-sorted CD34+ cells (>95% purity) and has GSM sample number GSM575. 

## Getting the data we want

Now that we have our GSE ID we can use the GEOquery to download the corresponding data as follows,


```r
# Download the mapping information and processed data This
# returns a list of eSets
GSE29617_set <- getGEO("GSE29617", destdir = "Data/GEO/")[[1]]
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29617/matrix/
## Found 1 file(s)
## GSE29617_series_matrix.txt.gz
## File stored at: 
## Data/GEO//GPL13158.soft
```
which returns (a list of) an ExpressionSet (eSet).

## The eSet class

What is an eSet? An S4 class that tries to:
- Coordinate high through-put (e.g., gene expression) and phenotype data.
- Provide common data container for diverse Bioconductor packages.


```r
str(GSE29617_set, max.level = 2)
```

```
## Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
##   ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
##   ..@ assayData        :<environment: 0x7ff331928540> 
##   ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
##   ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
##   ..@ annotation       : chr "GPL13158"
##   ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
##   ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
```

str() is the command to get the internal structure of an R object. 
An eSet contains the necessary "parts" to summarize an experiment.

## Classes and methods

**Everything in R is an OBJECT.**


- A class is the definition of an object.
- A method is a function that performs specific calculations on objects of a
specific class. Generic functions are used to determine the class of its
arguments and select the appropriate method. A generic function is a
function with a collection of methods.
- See ?Classes and ?Methods for more information.

## Classes and methods


```r
data(iris)
class(iris)
```

```
## [1] "data.frame"
```

```r
summary(iris)
```

```
##   Sepal.Length    Sepal.Width     Petal.Length    Petal.Width   
##  Min.   :4.300   Min.   :2.000   Min.   :1.000   Min.   :0.100  
##  1st Qu.:5.100   1st Qu.:2.800   1st Qu.:1.600   1st Qu.:0.300  
##  Median :5.800   Median :3.000   Median :4.350   Median :1.300  
##  Mean   :5.843   Mean   :3.057   Mean   :3.758   Mean   :1.199  
##  3rd Qu.:6.400   3rd Qu.:3.300   3rd Qu.:5.100   3rd Qu.:1.800  
##  Max.   :7.900   Max.   :4.400   Max.   :6.900   Max.   :2.500  
##        Species  
##  setosa    :50  
##  versicolor:50  
##  virginica :50  
##                 
##                 
## 
```

## Classes and methods

There are two types of classes in R: S3 Classes (old style, informal) and S4 Classes - (new style, more rigorous and formal)


```r
# S3 class
head(methods(class = "data.frame"))
```

```
## [1] "[.data.frame"    "[[.data.frame"   "[[<-.data.frame" "[<-.data.frame" 
## [5] "$.data.frame"    "$<-.data.frame"
```

```r
# S4 class
showMethods(class = "eSet")
```

```
## 
## Function ".DollarNames":
##  <not an S4 generic function>
## 
## Function ".SEW<-":
##  <not an S4 generic function>
## Function: [ (package base)
## x="eSet"
## 
## Function: [[ (package base)
## x="eSet"
## 
## Function: [[<- (package base)
## x="eSet"
## 
## 
## Function "%in%":
##  <not an S4 generic function>
## Function: $ (package base)
## x="eSet"
## 
## Function: $<- (package base)
## x="eSet"
## 
## Function: abstract (package Biobase)
## object="eSet"
## 
## 
## Function "active":
##  <not an S4 generic function>
## 
## Function "active<-":
##  <not an S4 generic function>
## 
## Function "addAttributes":
##  <not an S4 generic function>
## 
## Function "aggregate":
##  <not an S4 generic function>
## 
## Function "AIC":
##  <not an S4 generic function>
## Function: annotation (package BiocGenerics)
## object="eSet"
## 
## Function: annotation<- (package BiocGenerics)
## object="eSet", value="character"
## 
## 
## Function "applyFun":
##  <not an S4 generic function>
## 
## Function "applyFun<-":
##  <not an S4 generic function>
## 
## Function "applyParams":
##  <not an S4 generic function>
## 
## Function "applyParams<-":
##  <not an S4 generic function>
## 
## Function "as.env":
##  <not an S4 generic function>
## 
## Function "as.factor":
##  <not an S4 generic function>
## 
## Function "as.matrix":
##  <not an S4 generic function>
## 
## Function "as.table":
##  <not an S4 generic function>
## Function: assayData (package Biobase)
## object="eSet"
## 
## Function: assayData<- (package Biobase)
## object="eSet", value="AssayData"
## 
## 
## Function "BIC":
##  <not an S4 generic function>
## 
## Function "by":
##  <not an S4 generic function>
## 
## Function "chartr":
##  <not an S4 generic function>
## 
## Function "classNameForDisplay":
##  <not an S4 generic function>
## 
## Function "clearMemoryManagement":
##  <not an S4 generic function>
## 
## Function "clone":
##  <not an S4 generic function>
## 
## Function "coef":
##  <not an S4 generic function>
## Function: coerce (package methods)
## from="eSet", to="ExpressionSet"
## from="eSet", to="MultiSet"
## from="ExpressionSet", to="eSet"
## 
## 
## Function "collapse":
##  <not an S4 generic function>
## 
## Function "colnames<-":
##  <not an S4 generic function>
## 
## Function "columnMetadata":
##  <not an S4 generic function>
## 
## Function "columnMetadata<-":
##  <not an S4 generic function>
## Function: combine (package BiocGenerics)
## x="eSet", y="eSet"
## 
## 
## Function "comment.SAX":
##  <not an S4 generic function>
## 
## Function "compare":
##  <not an S4 generic function>
## 
## Function "compareRecursively":
##  <not an S4 generic function>
## 
## Function "complete.cases":
##  <not an S4 generic function>
## 
## Function "confint":
##  <not an S4 generic function>
## 
## Function "cor":
##  <not an S4 generic function>
## 
## Function "countMatches":
##  <not an S4 generic function>
## 
## Function "countOverlaps":
##  <not an S4 generic function>
## 
## Function "countQueryHits":
##  <not an S4 generic function>
## 
## Function "countSubjectHits":
##  <not an S4 generic function>
## 
## Function "cov":
##  <not an S4 generic function>
## 
## Function "coverage":
##  <not an S4 generic function>
## 
## Function "dbConn":
##  <not an S4 generic function>
## 
## Function "dbiDataType":
##  <not an S4 generic function>
## 
## Function "desc":
##  <not an S4 generic function>
## 
## Function "desc<-":
##  <not an S4 generic function>
## Function: description (package Biobase)
## object="eSet"
## 
## Function: description<- (package Biobase)
## object="eSet", value="MIAME"
## 
## 
## Function "diff":
##  <not an S4 generic function>
## Function: dim (package base)
## x="eSet"
## 
## Function: dimnames (package base)
## x="eSet"
## 
## Function: dimnames<- (package base)
## x="eSet"
## 
## Function: dims (package Biobase)
## object="eSet"
## 
## 
## Function "disjoin":
##  <not an S4 generic function>
## 
## Function "disjointBins":
##  <not an S4 generic function>
## 
## Function "distance":
##  <not an S4 generic function>
## 
## Function "distanceToNearest":
##  <not an S4 generic function>
## 
## Function "docName":
##  <not an S4 generic function>
## 
## Function "docName<-":
##  <not an S4 generic function>
## 
## Function "drop":
##  <not an S4 generic function>
## 
## Function "droplevels":
##  <not an S4 generic function>
## 
## Function "elementLengths":
##  <not an S4 generic function>
## 
## Function "elementMetadata":
##  <not an S4 generic function>
## 
## Function "elementMetadata<-":
##  <not an S4 generic function>
## 
## Function "elementType":
##  <not an S4 generic function>
## 
## Function "end":
##  <not an S4 generic function>
## 
## Function "end<-":
##  <not an S4 generic function>
## 
## Function "endElement.SAX":
##  <not an S4 generic function>
## 
## Function "endoapply":
##  <not an S4 generic function>
## 
## Function "entityDeclaration.SAX":
##  <not an S4 generic function>
## 
## Function "evalSeparately":
##  <not an S4 generic function>
## 
## Function "expand":
##  <not an S4 generic function>
## Function: experimentData (package Biobase)
## object="eSet"
## 
## Function: experimentData<- (package Biobase)
## object="eSet", value="MIAME"
## 
## 
## Function "extractList":
##  <not an S4 generic function>
## 
## Function "extractROWS":
##  <not an S4 generic function>
## Function: fData (package Biobase)
## object="eSet"
## 
## Function: fData<- (package Biobase)
## object="eSet", value="data.frame"
## 
## Function: featureData (package Biobase)
## object="eSet"
## 
## Function: featureData<- (package Biobase)
## object="eSet", value="AnnotatedDataFrame"
## 
## Function: featureNames (package Biobase)
## object="eSet"
## 
## Function: featureNames<- (package Biobase)
## object="eSet"
## 
## 
## Function "filterRules":
##  <not an S4 generic function>
## 
## Function "filterRules<-":
##  <not an S4 generic function>
## 
## Function "findMatches":
##  <not an S4 generic function>
## 
## Function "findOverlaps":
##  <not an S4 generic function>
## 
## Function "findRange":
##  <not an S4 generic function>
## 
## Function "findRun":
##  <not an S4 generic function>
## 
## Function "findXIncludeStartNodes":
##  <not an S4 generic function>
## 
## Function "fixedColumnNames":
##  <not an S4 generic function>
## 
## Function "flank":
##  <not an S4 generic function>
## 
## Function "flatten":
##  <not an S4 generic function>
## 
## Function "follow":
##  <not an S4 generic function>
## 
## Function "formals<-":
##  <not an S4 generic function>
## 
## Function "free":
##  <not an S4 generic function>
## 
## Function "functions":
##  <not an S4 generic function>
## Function: fvarLabels (package Biobase)
## object="eSet"
## 
## Function: fvarLabels<- (package Biobase)
## object="eSet"
## 
## Function: fvarMetadata (package Biobase)
## object="eSet"
## 
## Function: fvarMetadata<- (package Biobase)
## object="eSet", value="data.frame"
## 
## 
## Function "gaps":
##  <not an S4 generic function>
## 
## Function "getEffectiveNamespaces":
##  <not an S4 generic function>
## 
## Function "getEncoding":
##  <not an S4 generic function>
## 
## Function "getEncodingREnum":
##  <not an S4 generic function>
## 
## Function "getListElement":
##  <not an S4 generic function>
## 
## Function "getTableElementType":
##  <not an S4 generic function>
## 
## Function "GPL":
##  <not an S4 generic function>
## 
## Function "grouplength":
##  <not an S4 generic function>
## 
## Function "grouprank":
##  <not an S4 generic function>
## 
## Function "GSM":
##  <not an S4 generic function>
## 
## Function "gsub":
##  <not an S4 generic function>
## 
## Function "high2low":
##  <not an S4 generic function>
## 
## Function "ifelse":
##  <not an S4 generic function>
## Function: initialize (package methods)
## .Object="eSet"
## 
## 
## Function "IQR":
##  <not an S4 generic function>
## 
## Function "isConstant":
##  <not an S4 generic function>
## 
## Function "isDisjoint":
##  <not an S4 generic function>
## 
## Function "isEmpty":
##  <not an S4 generic function>
## 
## Function "isNormal":
##  <not an S4 generic function>
## 
## Function "iteratorFun":
##  <not an S4 generic function>
## 
## Function "iteratorFun<-":
##  <not an S4 generic function>
## 
## Function "levels":
##  <not an S4 generic function>
## 
## Function "Lfilter":
##  <not an S4 generic function>
## 
## Function "logLik":
##  <not an S4 generic function>
## 
## Function "low2high":
##  <not an S4 generic function>
## 
## Function "Ltablename":
##  <not an S4 generic function>
## 
## Function "mad":
##  <not an S4 generic function>
## 
## Function "map":
##  <not an S4 generic function>
## 
## Function "maskedratio":
##  <not an S4 generic function>
## 
## Function "maskedwidth":
##  <not an S4 generic function>
## 
## Function "mcols":
##  <not an S4 generic function>
## 
## Function "mcols<-":
##  <not an S4 generic function>
## 
## Function "mean":
##  <not an S4 generic function>
## 
## Function "median":
##  <not an S4 generic function>
## 
## Function "members":
##  <not an S4 generic function>
## 
## Function "mendoapply":
##  <not an S4 generic function>
## 
## Function "merge":
##  <not an S4 generic function>
## 
## Function "metadata<-":
##  <not an S4 generic function>
## 
## Function "mid":
##  <not an S4 generic function>
## 
## Function "mstack":
##  <not an S4 generic function>
## 
## Function "na.exclude":
##  <not an S4 generic function>
## 
## Function "na.omit":
##  <not an S4 generic function>
## 
## Function "narrow":
##  <not an S4 generic function>
## 
## Function "nchar":
##  <not an S4 generic function>
## 
## Function "nearest":
##  <not an S4 generic function>
## 
## Function "ngap":
##  <not an S4 generic function>
## 
## Function "nir_list":
##  <not an S4 generic function>
## 
## Function "nobj":
##  <not an S4 generic function>
## 
## Function "nobs":
##  <not an S4 generic function>
## 
## Function "normalizeSingleBracketReplacementValue":
##  <not an S4 generic function>
## Function: notes (package Biobase)
## object="eSet"
## 
## Function: notes<- (package Biobase)
## object="eSet", value="ANY"
## 
## 
## Function "nrun":
##  <not an S4 generic function>
## 
## Function "orderAsRanges":
##  <not an S4 generic function>
## 
## Function "overlapsAny":
##  <not an S4 generic function>
## 
## Function "packageName":
##  <not an S4 generic function>
## 
## Function "params":
##  <not an S4 generic function>
## Function: pData (package Biobase)
## object="eSet"
## 
## Function: pData<- (package Biobase)
## object="eSet", value="data.frame"
## 
## 
## Function "pgap":
##  <not an S4 generic function>
## Function: phenoData (package Biobase)
## object="eSet"
## 
## Function: phenoData<- (package Biobase)
## object="eSet", value="AnnotatedDataFrame"
## 
## 
## Function "pintersect":
##  <not an S4 generic function>
## 
## Function "plot":
##  <not an S4 generic function>
## 
## Function "pmap":
##  <not an S4 generic function>
## 
## Function "pop":
##  <not an S4 generic function>
## 
## Function "precede":
##  <not an S4 generic function>
## Function: preproc (package Biobase)
## object="eSet"
## 
## Function: preproc<- (package Biobase)
## object="eSet"
## 
## 
## Function "processingInstruction.SAX":
##  <not an S4 generic function>
## 
## Function "processSelfMatching":
##  <not an S4 generic function>
## 
## Function "profile":
##  <not an S4 generic function>
## 
## Function "promoters":
##  <not an S4 generic function>
## 
## Function "prompt":
##  <not an S4 generic function>
## Function: protocolData (package Biobase)
## object="eSet"
## 
## Function: protocolData<- (package Biobase)
## object="eSet", value="AnnotatedDataFrame"
## 
## 
## Function "psetdiff":
##  <not an S4 generic function>
## Function: pubMedIds (package Biobase)
## object="eSet"
## 
## Function: pubMedIds<- (package Biobase)
## object="eSet", value="character"
## 
## 
## Function "punion":
##  <not an S4 generic function>
## 
## Function "push":
##  <not an S4 generic function>
## 
## Function "quantile":
##  <not an S4 generic function>
## 
## Function "queryHits":
##  <not an S4 generic function>
## 
## Function "queryLength":
##  <not an S4 generic function>
## 
## Function "rangedData":
##  <not an S4 generic function>
## 
## Function "rangedData<-":
##  <not an S4 generic function>
## 
## Function "ranges":
##  <not an S4 generic function>
## 
## Function "ranges<-":
##  <not an S4 generic function>
## 
## Function "rdapply":
##  <not an S4 generic function>
## 
## Function "readHTMLList":
##  <not an S4 generic function>
## 
## Function "readHTMLTable":
##  <not an S4 generic function>
## 
## Function "readKeyValueDB":
##  <not an S4 generic function>
## 
## Function "readSolrDoc":
##  <not an S4 generic function>
## 
## Function "reduce":
##  <not an S4 generic function>
## 
## Function "reducerFun":
##  <not an S4 generic function>
## 
## Function "reducerFun<-":
##  <not an S4 generic function>
## 
## Function "reducerParams":
##  <not an S4 generic function>
## 
## Function "reducerParams<-":
##  <not an S4 generic function>
## 
## Function "reflect":
##  <not an S4 generic function>
## 
## Function "relistToClass":
##  <not an S4 generic function>
## 
## Function "removeAttributes":
##  <not an S4 generic function>
## 
## Function "removeXMLNamespaces":
##  <not an S4 generic function>
## 
## Function "rename":
##  <not an S4 generic function>
## 
## Function "replaceROWS":
##  <not an S4 generic function>
## 
## Function "reset":
##  <not an S4 generic function>
## 
## Function "resize":
##  <not an S4 generic function>
## 
## Function "restrict":
##  <not an S4 generic function>
## 
## Function "rev":
##  <not an S4 generic function>
## 
## Function "revElements":
##  <not an S4 generic function>
## 
## Function "reverse":
##  <not an S4 generic function>
## 
## Function "Rfilter":
##  <not an S4 generic function>
## 
## Function "Rle":
##  <not an S4 generic function>
## 
## Function "rownames<-":
##  <not an S4 generic function>
## 
## Function "Rtablename":
##  <not an S4 generic function>
## 
## Function "runLength":
##  <not an S4 generic function>
## 
## Function "runLength<-":
##  <not an S4 generic function>
## 
## Function "runmean":
##  <not an S4 generic function>
## 
## Function "runmed":
##  <not an S4 generic function>
## 
## Function "runq":
##  <not an S4 generic function>
## 
## Function "runsum":
##  <not an S4 generic function>
## 
## Function "runValue":
##  <not an S4 generic function>
## 
## Function "runValue<-":
##  <not an S4 generic function>
## 
## Function "runwtsum":
##  <not an S4 generic function>
## Function: sampleNames (package Biobase)
## object="eSet"
## 
## Function: sampleNames<- (package Biobase)
## object="eSet", value="ANY"
## 
## 
## Function "saveXML":
##  <not an S4 generic function>
## 
## Function "score":
##  <not an S4 generic function>
## 
## Function "score<-":
##  <not an S4 generic function>
## 
## Function "sd":
##  <not an S4 generic function>
## 
## Function "selectSomeIndex":
##  <not an S4 generic function>
## 
## Function "selfmatch":
##  <not an S4 generic function>
## 
## Function "setListElement":
##  <not an S4 generic function>
## 
## Function "shift":
##  <not an S4 generic function>
## 
## Function "shiftApply":
##  <not an S4 generic function>
## Function: show (package methods)
## object="eSet"
## 
## 
## Function "showAsCell":
##  <not an S4 generic function>
## 
## Function "simplify":
##  <not an S4 generic function>
## 
## Function "simplify<-":
##  <not an S4 generic function>
## 
## Function "simplifyNamespaces":
##  <not an S4 generic function>
## 
## Function "slice":
##  <not an S4 generic function>
## 
## Function "smoothEnds":
##  <not an S4 generic function>
## 
## Function "source":
##  <not an S4 generic function>
## 
## Function "space":
##  <not an S4 generic function>
## 
## Function "split":
##  <not an S4 generic function>
## 
## Function "split<-":
##  <not an S4 generic function>
## 
## Function "splitAsListReturnedClass":
##  <not an S4 generic function>
## 
## Function "splitRanges":
##  <not an S4 generic function>
## 
## Function "stack":
##  <not an S4 generic function>
## 
## Function "start":
##  <not an S4 generic function>
## 
## Function "start<-":
##  <not an S4 generic function>
## 
## Function "startElement.SAX":
##  <not an S4 generic function>
## Function: storageMode (package Biobase)
## object="eSet"
## 
## Function: storageMode<- (package Biobase)
## object="eSet", value="character"
## 
## 
## Function "sub":
##  <not an S4 generic function>
## 
## Function "subject":
##  <not an S4 generic function>
## 
## Function "subjectHits":
##  <not an S4 generic function>
## 
## Function "subjectLength":
##  <not an S4 generic function>
## 
## Function "subsetByFilter":
##  <not an S4 generic function>
## 
## Function "subsetByOverlaps":
##  <not an S4 generic function>
## 
## Function "substr":
##  <not an S4 generic function>
## 
## Function "substring":
##  <not an S4 generic function>
## 
## Function "subviews":
##  <not an S4 generic function>
## 
## Function "t":
##  <not an S4 generic function>
## 
## Function "text.SAX":
##  <not an S4 generic function>
## 
## Function "threebands":
##  <not an S4 generic function>
## 
## Function "tile":
##  <not an S4 generic function>
## 
## Function "togroup":
##  <not an S4 generic function>
## 
## Function "togrouplength":
##  <not an S4 generic function>
## 
## Function "togrouprank":
##  <not an S4 generic function>
## 
## Function "toHTML":
##  <not an S4 generic function>
## 
## Function "toList":
##  <not an S4 generic function>
## 
## Function "toLList":
##  <not an S4 generic function>
## 
## Function "tolower":
##  <not an S4 generic function>
## 
## Function "toRList":
##  <not an S4 generic function>
## 
## Function "toString":
##  <not an S4 generic function>
## 
## Function "toupper":
##  <not an S4 generic function>
## 
## Function "trim":
##  <not an S4 generic function>
## 
## Function "universe":
##  <not an S4 generic function>
## 
## Function "universe<-":
##  <not an S4 generic function>
## 
## Function "unsplit":
##  <not an S4 generic function>
## 
## Function "unstrsplit":
##  <not an S4 generic function>
## 
## Function "update":
##  <not an S4 generic function>
## Function: updateObject (package BiocGenerics)
## object="eSet"
## 
## Function: updateObjectTo (package Biobase)
## object="eSet", template="eSet"
## 
## 
## Function "values":
##  <not an S4 generic function>
## 
## Function "values<-":
##  <not an S4 generic function>
## 
## Function "var":
##  <not an S4 generic function>
## Function: varLabels (package Biobase)
## object="eSet"
## 
## Function: varLabels<- (package Biobase)
## object="eSet"
## 
## Function: varMetadata (package Biobase)
## object="eSet"
## 
## Function: varMetadata<- (package Biobase)
## object="eSet", value="data.frame"
## 
## 
## Function "vcov":
##  <not an S4 generic function>
## 
## Function "viewApply":
##  <not an S4 generic function>
## 
## Function "viewMaxs":
##  <not an S4 generic function>
## 
## Function "viewMeans":
##  <not an S4 generic function>
## 
## Function "viewMins":
##  <not an S4 generic function>
## 
## Function "viewRangeMaxs":
##  <not an S4 generic function>
## 
## Function "viewRangeMins":
##  <not an S4 generic function>
## 
## Function "Views":
##  <not an S4 generic function>
## 
## Function "viewSums":
##  <not an S4 generic function>
## 
## Function "viewWhichMaxs":
##  <not an S4 generic function>
## 
## Function "viewWhichMins":
##  <not an S4 generic function>
## 
## Function "vmembers":
##  <not an S4 generic function>
## 
## Function "which":
##  <not an S4 generic function>
## 
## Function "which.max":
##  <not an S4 generic function>
## 
## Function "which.min":
##  <not an S4 generic function>
## 
## Function "whichFirstNotNormal":
##  <not an S4 generic function>
## 
## Function "width":
##  <not an S4 generic function>
## 
## Function "width<-":
##  <not an S4 generic function>
## 
## Function "window":
##  <not an S4 generic function>
## 
## Function "window<-":
##  <not an S4 generic function>
## 
## Function "with":
##  <not an S4 generic function>
## 
## Function "within":
##  <not an S4 generic function>
## 
## Function "xmlAttrs<-":
##  <not an S4 generic function>
## 
## Function "xmlAttrsToDataFrame":
##  <not an S4 generic function>
## 
## Function "xmlChildren<-":
##  <not an S4 generic function>
## 
## Function "xmlClone":
##  <not an S4 generic function>
## 
## Function "xmlNamespace<-":
##  <not an S4 generic function>
## 
## Function "xmlNamespaces<-":
##  <not an S4 generic function>
## 
## Function "xmlParent":
##  <not an S4 generic function>
## 
## Function "xmlRoot<-":
##  <not an S4 generic function>
## 
## Function "xmlSource":
##  <not an S4 generic function>
## 
## Function "xmlSourceFunctions":
##  <not an S4 generic function>
## 
## Function "xmlSourceSection":
##  <not an S4 generic function>
## 
## Function "xmlSourceTask":
##  <not an S4 generic function>
## 
## Function "xmlSourceThread":
##  <not an S4 generic function>
## 
## Function "xmlToDataFrame":
##  <not an S4 generic function>
## 
## Function "xmlToS4":
##  <not an S4 generic function>
## 
## Function "xmlValue<-":
##  <not an S4 generic function>
```

## The eSet

You can get a sense of the defined methods for an `eSet` as follows:

```r
library(Biobase)
showMethods(class = "eSet")
```
in particular, the following methods are rather convenient:

- assayData(obj); assayData(obj) `<-` value: access or assign assayData
- phenoData(obj); phenoData(obj) `<-` value: access or assign phenoData
- experimentData(obj); experimentData(obj) `<-` value: access or assign experimentData
- annotation(obj); annotation(obj) `<-` value: access or assign annotation

## The ExpressionSet subclass

Similar to the `eSet` class but tailored to gene expression, with an expression matrix that can be accessed with the `exprs` method.


```r
class(GSE29617_set)
```

```
## [1] "ExpressionSet"
## attr(,"package")
## [1] "Biobase"
```

```r
exprs(GSE29617_set)[1:2, 1:3]
```

```
##              GSM733942 GSM733943 GSM733944
## 1007_PM_s_at   4.62965   4.46083   4.54702
## 1053_PM_at     4.41375   4.72417   4.48257
```

also provides additional methods such as `fData`.

## The ExpressionSet subclass

`ExpressionSet` objects are meant to facilitate the adoption of MIAME standard. MIAME = "Minimum Information about a Microarray experiment". Alvis Brazma et. al. (2001) Nature Genetics
Unfortrunately, not all contributors will upload all the information.

```r
# Information about preprocessing Nothing in here!
preproc(GSE29617_set)
```

```
## list()
```

## The ExpressionSet subclass


```r
# A data.frame with number of rows equal to the number of
# samples
pData(GSE29617_set)[1:2, 1:2]
```

```
##                                                  title geo_accession
## GSM733942 2008 TIV subject ID 2 at D0 post-vaccination     GSM733942
## GSM733943 2008 TIV subject ID 2 at D7 post-vaccination     GSM733943
```

```r
# A data.frame with number of rows equal to the number of
# features/probes
fData(GSE29617_set)[1:2, 1:2]
```

```
##                        ID GB_ACC
## 1007_PM_s_at 1007_PM_s_at U48705
## 1053_PM_at     1053_PM_at M87338
```

## The ExpressionSet subclass 

So the `ExpressionSet` objects facilitate the encapsulation of everything that's needed to summarize and analyze an experiment. Specific elements can be access with the `@` operator but many classes have convenient accessor methods.


```r
fData(GSE29617_set)[1:2, 1:2]
```

```
##                        ID GB_ACC
## 1007_PM_s_at 1007_PM_s_at U48705
## 1053_PM_at     1053_PM_at M87338
```

```r
# Note that S4 classes can be nested!
GSE29617_set@featureData@data[1:2, 1:2]
```

```
##                        ID GB_ACC
## 1007_PM_s_at 1007_PM_s_at U48705
## 1053_PM_at     1053_PM_at M87338
```

## What if you want the raw data?

GEO also provides access to raw data that can be downloaded with `GEOquery`.



```r
# Download all raw data. This should only be evaluated once
# Then the data would be stored locally in the data directory
# Make sure the directory exists
if (length(dir("Data/GEO/", pattern = "GSE29617")) == 0) {
    getGEOSuppFiles("GSE29617", makeDirectory = TRUE, baseDir = "./Data/GEO/")
    untar("./Data/GEO/GSE29617/GSE29617_RAW.tar", exdir = "./Data/GEO/GSE29617/", 
        tar = Sys.getenv("TAR"))
}
# untar downloaded data
```

## Starting from the raw data

Now that we have the Affymetrix raw data (CEL) files, we can apply some of the concepts we've discussed related to normalization and probe summary. We first need to load the appropriate package



```r
biocLite("affy")
```

```r
library(affy)
```

then we use the following commands

```r
# Read the CEL file and creates and AffyBatch
GSE29617_affyBatch <- ReadAffy(celfile.path = "Data/GEO/GSE29617/")
# Normalize and summarize the data
GSE29617_set2 <- rma(GSE29617_affyBatch)
```

```
## Background correcting
## Normalizing
## Calculating Expression
```

## Starting from the raw data

Let's check the results and compare to the expression matrix that was submitted to GEO

```r
exprs(GSE29617_set2)[1:2, 1:2]
```

```
##              GSM733942.CEL.gz GSM733943.CEL.gz
## 1007_PM_s_at         4.629650         4.460843
## 1053_PM_at           4.413773         4.724213
```

The rows are the features (i.e., probes). Columns are the samples.

## What are those probes?


```r
# We first need to install our annotation package
library(BiocInstaller)
# Note that you don't have to use source anymore!
biocLite("hthgu133a.db")
```



```r
library(hthgu133a.db)
probe_ids <- rownames(GSE29617_set2)
probe_data <- select(hthgu133a.db, keys = probe_ids, columns = "SYMBOL", 
    keytype = "PROBEID")
probe_data[1, ]
```

```
##        PROBEID SYMBOL
## 1 1007_PM_s_at   <NA>
```
This didn't work very well, did it?
The problem is that the probe IDs in hthgu133a.db have a different naming scheme than those in GSE29617_set2. This is fixed on the next slide.

## What are those probes?

Let's fix this: Replace _PM with <empty> for the probe id names in GSE29617_set2

```r
probe_ids <- gsub("_PM", "", rownames(GSE29617_set2))
probe_data <- select(hthgu133a.db, keys = probe_ids, columns = "SYMBOL", 
    keytype = "PROBEID")
```

```
## Warning in .generateExtraRows(tab, keys, jointype): 'select' resulted in
## 1:many mapping between keys and return rows
```

```r
probe_data[1, ]
```

```
##     PROBEID SYMBOL
## 1 1007_s_at   DDR1
```
What's the warning? Some probes match up with multiple genes, therefore those probe IDs will have more than one record.

## What are those probes?


This gives us too many rows, what do we do? Concatenate the gene names so that there will be one row per probe ID.


```r
library(data.table)
probe_data_dt <- data.table(probe_data)
probe_data_dt_unique <- probe_data_dt[, list(SYMBOL = paste(SYMBOL, 
    collapse = ";")), by = "PROBEID"]
probe_data_dt_unique[SYMBOL %like% ";"]
```

```
##           PROBEID                                  SYMBOL
##    1:   1007_s_at                            DDR1;MIR4640
##    2:     1294_at                            UBA7;MIR5193
##    3:     1773_at                        FNTB;CHURC1-FNTB
##    4: 200012_x_at         RPL21;SNORD102;SNORA27;RPL21P28
##    5:   200018_at             RPS13;SNORD14B;LOC100508408
##   ---                                                    
## 1143:  65133_i_at                      INO80B;INO80B-WBP1
## 1144:    65585_at FAM86C1;FAM86B1;FAM86FP;FAM86B2;FAM86DP
## 1145:    66053_at                 HNRNPUL2;HNRNPUL2-BSCL2
## 1146:    78495_at                        LOC155060;ZNF783
## 1147:    91617_at                           DGCR8;MIR1306
```

## Completing our ExpressionSet


```r
annotaded_probes <- data.frame(probe_data_dt_unique)
rownames(annotaded_probes) <- rownames(GSE29617_set2)
fData(GSE29617_set2) <- annotaded_probes
head(fData(GSE29617_set2))
```

```
##                PROBEID       SYMBOL
## 1007_PM_s_at 1007_s_at DDR1;MIR4640
## 1053_PM_at     1053_at         RFC2
## 117_PM_at       117_at        HSPA6
## 121_PM_at       121_at         PAX8
## 1255_PM_g_at 1255_g_at       GUCA1A
## 1294_PM_at     1294_at UBA7;MIR5193
```

Now you're ready to get the data that you will have to analyze for your second homework. 

