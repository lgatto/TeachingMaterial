# Indentification data using `mzR`, `mzID`, `MSnID`

## Handling identification data

There are two packages that can be used to parse `mzIdentML` files,
namely `mzR` (that we have already used for raw data) and `mzID`. The
major difference is that the former leverages C++ code from
`proteowizard` and is hence faster than the latter (which uses the
`XML` R package). They both work in similar ways.


|   |Data type      |File format |Data structure |Package |
|:--|:--------------|:-----------|:--------------|:-------|
|4  |Identification |mzIdentML   |mzRident       |mzR     |
|5  |Identification |mzIdentML   |mzID           |mzID    |

We are going to use the following identification file in this practical:


```r
library("msdata")
idf <- ident(full.names = TRUE)
idf
```

```
## [1] "/home/lg390/R/x86_64-pc-linux-gnu-library/3.3/msdata/ident/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzid"
```

### mzID

The main functions are `mzID` to read the data into a dedicated data
class and `flatten` to transform it into a `data.frame`. 


```r
library("mzID")
id <- mzID(idf)
```

```
## reading TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzid... DONE!
```

```r
id
```

```
## An mzID object
## 
## Software used:   MS-GF+ (version: Beta (v10072))
## 
## Rawfile:         /home/lg390/dev/01_svn/workflows/proteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML
## 
## Database:        /home/lg390/dev/01_svn/workflows/proteomics/erwinia_carotovora.fasta
## 
## Number of scans: 5343
## Number of PSM's: 5656
```

Various data can be extracted from the `mzID` object, using one the
accessor functions such as `database`, `sofware`, `scans`, `peptides`,
... The object can also be converted into a `data.frame` using the
`flatten` function.


```r
head(flatten(id))
```

```
## Error in head(flatten(id)): could not find function "flatten"
```

#### Exercise

Open the /home/lg390/R/x86_64-pc-linux-gnu-library/3.3/msdata/ident/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzid file as shown above. What software and how many
proteins were in the database used to perform the search?

### `mzR`

The `mzR` interface provides a similar interface. It is however much
faster as it does not read all the data into memory and only extracts
relevant data on demand. It has also accessor functions such as
`softwareInfo`, `mzidInfo`, ... (use `showMethods(classes = "mzRident", where = "package:mzR")`)
to see all available methods.


```r
library("mzR")
```

```
## Loading required package: Rcpp
```

```r
id2 <- openIDfile(idf)
id2
```

```
## Identification file handle.
## Filename:  TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzid 
## Number of psms:  5759
```

```r
softwareInfo(id2)
```

```
## [1] "MS-GF+ Beta (v10072) "                       
## [2] "ProteoWizard MzIdentML 3.0.9490 ProteoWizard"
```

The identification data can be accessed as a `data.frame` with the
`psms` accessor.


```r
head(psms(id2))
```

```
##                                      spectrumID chargeState rank
## 1 controllerType=0 controllerNumber=1 scan=5782           3    1
## 2 controllerType=0 controllerNumber=1 scan=6037           3    1
## 3 controllerType=0 controllerNumber=1 scan=5235           3    1
## 4 controllerType=0 controllerNumber=1 scan=5397           3    1
## 5 controllerType=0 controllerNumber=1 scan=6075           3    1
##   passThreshold experimentalMassToCharge calculatedMassToCharge
## 1          TRUE                1080.2325              1080.2321
## 2          TRUE                1002.2089              1002.2115
## 3          TRUE                1189.2836              1189.2800
## 4          TRUE                 960.5365               960.5365
## 5          TRUE                1264.3409              1264.3419
##                              sequence modNum isDecoy post pre start end
## 1 PVQIQAGEDSNVIGALGGAVLGGFLGNTIGGGSGR      0   FALSE    S   R    50  84
## 2        TQVLDGLINANDIEVPVALIDGEIDVLR      0   FALSE    R   K   288 315
## 3   TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR      0   FALSE    A   R   192 224
## 4         SQILQQAGTSVLSQANQVPQTVLSLLR      0   FALSE    -   R   264 290
## 5 PIIGDNPFVVVLPDVVLDESTADQTQENLALLISR      0   FALSE    F   R   119 153
##   DatabaseAccess DBseqLength DatabaseSeq
## 1        ECA1932         155            
## 2        ECA1147         434            
## 3        ECA0013         295            
## 4        ECA1731         290            
## 5        ECA1443         298            
##                                         DatabaseDescription acquisitionNum
## 1                        ECA1932 outer membrane lipoprotein           5782
## 2                                    ECA1147 trigger factor           6037
## 3                ECA0013 ribose-binding periplasmic protein           5235
## 4                                         ECA1731 flagellin           5397
## 5      ECA1443 UTP--glucose-1-phosphate uridylyltransferase           6075
##  [ reached getOption("max.print") -- omitted 1 row ]
```

#### Exercise

Is there a relation between the length of a protein and the number of
identified peptides, conditioned by the (average) e-value of the
identifications?



## Exploration and Assessment of Confidence of LC-MSn Proteomics Identifications using `MSnID`

Extracts MS/MS ID data from mzIdentML (leveraging the `mzID` package)
or text files. After collating the search results from multiple
datasets it assesses their identification quality and optimize
filtering criteria to achieve the maximum number of identifications
while not exceeding a specified false discovery rate. Also contains a
number of utilities to explore the MS/MS results and assess missed and
irregular enzymatic cleavages, mass measurement accuracy, etc.

Let's reproduce the analysis described the `MSnID` vignette. Open it
with


```r
vignette("msnid_vignette", package = "MSnID")
```
