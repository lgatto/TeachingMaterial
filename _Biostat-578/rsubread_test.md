# RNA-seq Pipeline: Rsubread, limma, and edgeR

## Case study: using a Bioconductor R pipeline to analyze RNA-seq data

This example has been adapted from the work of:

Wei Shi (shi at wehi dot edu dot au), Yang Liao and Gordon K Smyth
Bioinformatics Division, Walter and Eliza Hall Institute, Melbourne, Australia

- Case: http://bioinf.wehi.edu.au/RNAseqCaseStudy/
- Code: http://bioinf.wehi.edu.au/RNAseqCaseStudy/code.R
- Data: http://bioinf.wehi.edu.au/RNAseqCaseStudy/data.tar.gz

Requirements:

- The version of Rsubread package should be 1.12.1 or later. 
- You should run R version 3.0.2 or later. 

## Record the start time


```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:12:48 PST"
```

## Setup

> Load libraries, read in names of FASTQ files and create a design matrix.

### Set options

Set global knitr options.


```r
library("knitr")
opts_chunk$set(tidy=FALSE, cache=FALSE, messages=FALSE)
```

### Install packages


```r
packages <- c("Rsubread", "limma", "edgeR")
source("http://bioconductor.org/biocLite.R")
```

```
## Bioconductor version 3.0 (BiocInstaller 1.16.1), ?biocLite for help
```

```r
for (pkg in packages)
{
    require(pkg, character.only = TRUE, quietly = TRUE) || biocLite(pkg) 
}
```

### Load libraries


```r
library(Rsubread)
library(limma)
library(edgeR)
print(Sys.time())
```

```
## [1] "2015-03-05 12:12:50 PST"
```

### Download data file


```r
# Create the data folder if it does not already exist
datadir <- "./Data/rsubread_test"
dir.create(file.path(datadir), showWarnings=FALSE, recursive=TRUE)

# Enter data folder, first saving location of current folder
projdir <- getwd()
setwd(datadir)
datadir <- getwd()

# Get the file
dataurl <- "http://bioinf.wehi.edu.au/RNAseqCaseStudy/data.tar.gz"
datafile <- "data.tar.gz"
if (! file.exists(datafile)) {
    print("Downloading data file...")
    download.file(dataurl, datafile) # 282.3 Mb
}
```

```
## [1] "Downloading data file..."
```

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:19:17 PST"
```

### Extract data file


```r
setwd(datadir)
targetsfile <- "Targets.txt"
if (! file.exists(targetsfile)) {
    print("Extracting data file...")
    untar(datafile, tar="internal")
}
```

```
## [1] "Extracting data file..."
```

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:19:55 PST"
```

### Read target file


```r
setwd(datadir)
options(digits=2)
targets <- readTargets(targetsfile)
targets
```

```
##   CellType  InputFile OutputFile
## 1        A A_1.txt.gz    A_1.bam
## 2        A A_2.txt.gz    A_2.bam
## 3        B B_1.txt.gz    B_1.bam
## 4        B B_2.txt.gz    B_2.bam
```

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:19:55 PST"
```

### Create design matrix


```r
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)
print(Sys.time())
```

```
## [1] "2015-03-05 12:19:55 PST"
```

## Build reference

> Build an index for human chromosome 1. This typically takes ~3 minutes. 
Index files with basename 'chr1' will be generated in your current working 
directory. 


```r
setwd(datadir)
buildindex(basename="chr1",reference="hg19_chr1.fa")
```

```
## 
##         ==========     _____ _    _ ____  _____  ______          _____  
##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
##        Rsubread 1.16.1
## 
## //=========================== indexBuilder setting ===========================\\
## ||                                                                            ||
## ||                Index name : chr1                                           ||
## ||               Index space : base-space                                     ||
## ||                    Memory : 8000 Mbytes                                    ||
## ||          Repeat threshold : 24 repeats                                     ||
## ||  Distance to next subread : 3                                              ||
## ||                                                                            ||
## ||               Input files : 1 file in total                                ||
## ||                             o hg19_chr1.fa                                 ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## //================================= Running ==================================\\
## ||                                                                            ||
## || Check the integrity of provided reference sequences ...                    ||
## || No format issues were found                                                ||
## || Scan uninformative subreads in reference sequences ...                     ||
## ||    8%,   0 mins elapsed, rate=3532.5k bps/s, total=249m                    ||
## ||   16%,   0 mins elapsed, rate=3914.6k bps/s, total=249m                    ||
## ||   24%,   0 mins elapsed, rate=4087.1k bps/s, total=249m                    ||
## ||   33%,   0 mins elapsed, rate=4194.7k bps/s, total=249m                    ||
## ||   41%,   0 mins elapsed, rate=4247.1k bps/s, total=249m                    ||
## ||   49%,   0 mins elapsed, rate=4357.1k bps/s, total=249m                    ||
## ||   58%,   0 mins elapsed, rate=4787.0k bps/s, total=249m                    ||
## ||   66%,   0 mins elapsed, rate=4725.6k bps/s, total=249m                    ||
## ||   74%,   0 mins elapsed, rate=4670.6k bps/s, total=249m                    ||
## ||   83%,   0 mins elapsed, rate=4642.6k bps/s, total=249m                    ||
## ||   91%,   0 mins elapsed, rate=4606.1k bps/s, total=249m                    ||
## ||   99%,   1 mins elapsed, rate=4571.4k bps/s, total=249m                    ||
## || 64353 uninformative subreads were found.                                   ||
## || These subreads were excluded from index building.                          ||
## || Build the index...                                                         ||
## ||    8%,   1 mins elapsed, rate=1611.5k bps/s, total=249m                    ||
## ||   16%,   1 mins elapsed, rate=1779.0k bps/s, total=249m                    ||
## ||   24%,   1 mins elapsed, rate=1822.1k bps/s, total=249m                    ||
## ||   33%,   1 mins elapsed, rate=1835.8k bps/s, total=249m                    ||
## ||   41%,   1 mins elapsed, rate=1845.1k bps/s, total=249m                    ||
## ||   49%,   2 mins elapsed, rate=1894.3k bps/s, total=249m                    ||
## ||   58%,   2 mins elapsed, rate=2121.3k bps/s, total=249m                    ||
## ||   66%,   2 mins elapsed, rate=2101.3k bps/s, total=249m                    ||
## ||   74%,   2 mins elapsed, rate=2075.7k bps/s, total=249m                    ||
## ||   83%,   2 mins elapsed, rate=2056.4k bps/s, total=249m                    ||
## ||   91%,   2 mins elapsed, rate=2039.2k bps/s, total=249m                    ||
## ||   99%,   3 mins elapsed, rate=2027.6k bps/s, total=249m                    ||
## || Save current index block...                                                ||
## ||  [ 0.0% finished ]                                                         ||
## ||  [ 10.0% finished ]                                                        ||
## ||  [ 20.0% finished ]                                                        ||
## ||  [ 30.0% finished ]                                                        ||
## ||  [ 40.0% finished ]                                                        ||
## ||  [ 50.0% finished ]                                                        ||
## ||  [ 60.0% finished ]                                                        ||
## ||  [ 70.0% finished ]                                                        ||
## ||  [ 80.0% finished ]                                                        ||
## ||  [ 90.0% finished ]                                                        ||
## ||  [ 100.0% finished ]                                                       ||
## ||                                                                            ||
## ||                      Total running time: 3.3 minutes.                      ||
## ||                     Index chr1 was successfully built!                     ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
```

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:23:11 PST"
```

## Align reads

> Perform read alignment for all four libraries and report uniquely mapped 
reads only. This typically takes ~5 minutes. The generated SAM files, which 
include the mapping results, will be saved in your current working directory. 


```r
setwd(datadir)
align(index="chr1", readfile1=targets$InputFile, 
      input_format="gzFASTQ", output_format="BAM", 
      output_file=targets$OutputFile, 
      tieBreakHamming=TRUE, unique=TRUE, indels=5)
```

```
## 
##         ==========     _____ _    _ ____  _____  ______          _____  
##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
##        Rsubread 1.16.1
## 
## //========================== subread-align setting ===========================\\
## ||                                                                            ||
## ||           Function : Read alignment                                        ||
## ||            Threads : 1                                                     ||
## ||         Input file : A_1.txt.gz                                            ||
## ||        Output file : A_1.bam (BAM)                                         ||
## ||         Index name : chr1                                                  ||
## ||       Phred offset : 33                                                    ||
## ||                                                                            ||
## ||          Min votes : 3                                                     ||
## ||         Max indels : 5                                                     ||
## ||  # of Best mapping : 1                                                     ||
## ||     Unique mapping : yes                                                   ||
## ||   Hamming distance : yes                                                   ||
## ||     Quality scores : no                                                    ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## //====================== Running (05-Mar-2015 12:23:11) ======================\\
## ||                                                                            ||
## || Decompress A_1.txt.gz...                                                   ||
## || The input file contains base space reads.                                  ||
## || Load the 1-th index block...                                               ||
## || Map reads...                                                               ||
## ||    0% completed,   0 mins elapsed, total=657k reads, rate=9.2k/s           ||
## ||    6% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   12% completed,   0 mins elapsed, total=657k reads, rate=9.2k/s           ||
## ||   18% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   24% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   30% completed,   0 mins elapsed, total=657k reads, rate=9.2k/s           ||
## ||   36% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   42% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   48% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   54% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## || Detect indels...                                                           ||
## ||   60% completed,   0 mins elapsed, total=657k reads, rate=8.6k/s           ||
## ||   62% completed,   0 mins elapsed, total=657k reads, rate=8.7k/s           ||
## ||   63% completed,   0 mins elapsed, total=657k reads, rate=8.8k/s           ||
## ||   64% completed,   0 mins elapsed, total=657k reads, rate=8.9k/s           ||
## ||   65% completed,   0 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   66% completed,   0 mins elapsed, total=657k reads, rate=9.1k/s           ||
## ||   67% completed,   0 mins elapsed, total=657k reads, rate=9.2k/s           ||
## ||   68% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   69% completed,   0 mins elapsed, total=657k reads, rate=9.4k/s           ||
## || Realign reads...                                                           ||
## ||   71% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   72% completed,   0 mins elapsed, total=657k reads, rate=9.3k/s           ||
## ||   73% completed,   0 mins elapsed, total=657k reads, rate=9.2k/s           ||
## ||   74% completed,   0 mins elapsed, total=657k reads, rate=9.2k/s           ||
## ||   75% completed,   0 mins elapsed, total=657k reads, rate=9.1k/s           ||
## ||   76% completed,   0 mins elapsed, total=657k reads, rate=9.1k/s           ||
## ||   77% completed,   0 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   78% completed,   0 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   79% completed,   0 mins elapsed, total=657k reads, rate=9.0k/s           ||
## || 657474 reads were processed. Save the mapping results for them...          ||
## ||   81% completed,   1 mins elapsed, total=657k reads, rate=8.9k/s           ||
## ||   84% completed,   1 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   86% completed,   1 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   88% completed,   1 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   90% completed,   1 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   92% completed,   1 mins elapsed, total=657k reads, rate=9.0k/s           ||
## ||   94% completed,   1 mins elapsed, total=657k reads, rate=9.1k/s           ||
## ||   96% completed,   1 mins elapsed, total=657k reads, rate=9.1k/s           ||
## ||   98% completed,   1 mins elapsed, total=657k reads, rate=9.1k/s           ||
## ||                                                                            ||
## ||                          Completed successfully.                           ||
## ||                                                                            ||
## \\============================================================================//
## 
## //================================= Summary ==================================\\
## ||                                                                            ||
## ||          Processed : 657474 reads                                          ||
## ||             Mapped : 609482 reads (92.7%)                                  ||
## ||             Indels : 4099                                                  ||
## ||                                                                            ||
## ||       Running time : 1.2 minutes                                           ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## 
##         ==========     _____ _    _ ____  _____  ______          _____  
##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
##        Rsubread 1.16.1
## 
## //========================== subread-align setting ===========================\\
## ||                                                                            ||
## ||           Function : Read alignment                                        ||
## ||            Threads : 1                                                     ||
## ||         Input file : A_2.txt.gz                                            ||
## ||        Output file : A_2.bam (BAM)                                         ||
## ||         Index name : chr1                                                  ||
## ||       Phred offset : 33                                                    ||
## ||                                                                            ||
## ||          Min votes : 3                                                     ||
## ||         Max indels : 5                                                     ||
## ||  # of Best mapping : 1                                                     ||
## ||     Unique mapping : yes                                                   ||
## ||   Hamming distance : yes                                                   ||
## ||     Quality scores : no                                                    ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## //====================== Running (05-Mar-2015 12:24:23) ======================\\
## ||                                                                            ||
## || Decompress A_2.txt.gz...                                                   ||
## || The input file contains base space reads.                                  ||
## || Load the 1-th index block...                                               ||
## || Map reads...                                                               ||
## ||    1% completed,   0 mins elapsed, total=562k reads, rate=9.7k/s           ||
## ||    7% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   13% completed,   0 mins elapsed, total=562k reads, rate=9.8k/s           ||
## ||   19% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   25% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   31% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   37% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   43% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   49% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   55% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## || Detect indels...                                                           ||
## ||   61% completed,   0 mins elapsed, total=562k reads, rate=9.2k/s           ||
## ||   62% completed,   0 mins elapsed, total=562k reads, rate=9.3k/s           ||
## ||   63% completed,   0 mins elapsed, total=562k reads, rate=9.4k/s           ||
## ||   64% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   65% completed,   0 mins elapsed, total=562k reads, rate=9.6k/s           ||
## ||   66% completed,   0 mins elapsed, total=562k reads, rate=9.7k/s           ||
## ||   67% completed,   0 mins elapsed, total=562k reads, rate=9.8k/s           ||
## ||   68% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   69% completed,   0 mins elapsed, total=562k reads, rate=10.0k/s          ||
## || Realign reads...                                                           ||
## ||   71% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   72% completed,   0 mins elapsed, total=562k reads, rate=9.9k/s           ||
## ||   73% completed,   0 mins elapsed, total=562k reads, rate=9.8k/s           ||
## ||   74% completed,   0 mins elapsed, total=562k reads, rate=9.8k/s           ||
## ||   75% completed,   0 mins elapsed, total=562k reads, rate=9.7k/s           ||
## ||   76% completed,   0 mins elapsed, total=562k reads, rate=9.7k/s           ||
## ||   77% completed,   0 mins elapsed, total=562k reads, rate=9.6k/s           ||
## ||   78% completed,   0 mins elapsed, total=562k reads, rate=9.6k/s           ||
## ||   79% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## || 562490 reads were processed. Save the mapping results for them...          ||
## ||   82% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   84% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   86% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   88% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   90% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   92% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   94% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   96% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||   98% completed,   0 mins elapsed, total=562k reads, rate=9.5k/s           ||
## ||                                                                            ||
## ||                          Completed successfully.                           ||
## ||                                                                            ||
## \\============================================================================//
## 
## //================================= Summary ==================================\\
## ||                                                                            ||
## ||          Processed : 562490 reads                                          ||
## ||             Mapped : 521918 reads (92.8%)                                  ||
## ||             Indels : 3511                                                  ||
## ||                                                                            ||
## ||       Running time : 1.0 minutes                                           ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## 
##         ==========     _____ _    _ ____  _____  ______          _____  
##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
##        Rsubread 1.16.1
## 
## //========================== subread-align setting ===========================\\
## ||                                                                            ||
## ||           Function : Read alignment                                        ||
## ||            Threads : 1                                                     ||
## ||         Input file : B_1.txt.gz                                            ||
## ||        Output file : B_1.bam (BAM)                                         ||
## ||         Index name : chr1                                                  ||
## ||       Phred offset : 33                                                    ||
## ||                                                                            ||
## ||          Min votes : 3                                                     ||
## ||         Max indels : 5                                                     ||
## ||  # of Best mapping : 1                                                     ||
## ||     Unique mapping : yes                                                   ||
## ||   Hamming distance : yes                                                   ||
## ||     Quality scores : no                                                    ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## //====================== Running (05-Mar-2015 12:25:22) ======================\\
## ||                                                                            ||
## || Decompress B_1.txt.gz...                                                   ||
## || The input file contains base space reads.                                  ||
## || Load the 1-th index block...                                               ||
## || Map reads...                                                               ||
## ||    0% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||    6% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   12% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   18% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   24% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   30% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   36% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   42% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   48% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   54% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## || Detect indels...                                                           ||
## ||   61% completed,   0 mins elapsed, total=902k reads, rate=9.6k/s           ||
## ||   62% completed,   0 mins elapsed, total=902k reads, rate=9.7k/s           ||
## ||   63% completed,   0 mins elapsed, total=902k reads, rate=9.8k/s           ||
## ||   64% completed,   0 mins elapsed, total=902k reads, rate=9.9k/s           ||
## ||   65% completed,   0 mins elapsed, total=902k reads, rate=10.0k/s          ||
## ||   66% completed,   0 mins elapsed, total=902k reads, rate=10.1k/s          ||
## ||   67% completed,   0 mins elapsed, total=902k reads, rate=10.2k/s          ||
## ||   68% completed,   0 mins elapsed, total=902k reads, rate=10.3k/s          ||
## ||   69% completed,   0 mins elapsed, total=902k reads, rate=10.4k/s          ||
## || Realign reads...                                                           ||
## ||   71% completed,   1 mins elapsed, total=902k reads, rate=10.1k/s          ||
## ||   72% completed,   1 mins elapsed, total=902k reads, rate=9.9k/s           ||
## ||   73% completed,   1 mins elapsed, total=902k reads, rate=9.6k/s           ||
## ||   74% completed,   1 mins elapsed, total=902k reads, rate=9.4k/s           ||
## ||   75% completed,   1 mins elapsed, total=902k reads, rate=9.2k/s           ||
## ||   76% completed,   1 mins elapsed, total=902k reads, rate=9.1k/s           ||
## ||   77% completed,   1 mins elapsed, total=902k reads, rate=8.9k/s           ||
## ||   78% completed,   1 mins elapsed, total=902k reads, rate=8.7k/s           ||
## ||   79% completed,   1 mins elapsed, total=902k reads, rate=8.6k/s           ||
## || 902841 reads were processed. Save the mapping results for them...          ||
## ||   82% completed,   1 mins elapsed, total=902k reads, rate=8.5k/s           ||
## ||   84% completed,   1 mins elapsed, total=902k reads, rate=8.5k/s           ||
## ||   86% completed,   1 mins elapsed, total=902k reads, rate=8.5k/s           ||
## ||   88% completed,   1 mins elapsed, total=902k reads, rate=8.6k/s           ||
## ||   90% completed,   1 mins elapsed, total=902k reads, rate=8.6k/s           ||
## ||   92% completed,   1 mins elapsed, total=902k reads, rate=8.6k/s           ||
## ||   94% completed,   1 mins elapsed, total=902k reads, rate=8.7k/s           ||
## ||   96% completed,   1 mins elapsed, total=902k reads, rate=8.7k/s           ||
## ||   98% completed,   1 mins elapsed, total=902k reads, rate=8.7k/s           ||
## ||                                                                            ||
## ||                          Completed successfully.                           ||
## ||                                                                            ||
## \\============================================================================//
## 
## //================================= Summary ==================================\\
## ||                                                                            ||
## ||          Processed : 902841 reads                                          ||
## ||             Mapped : 838835 reads (92.9%)                                  ||
## ||             Indels : 5024                                                  ||
## ||                                                                            ||
## ||       Running time : 1.7 minutes                                           ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## 
##         ==========     _____ _    _ ____  _____  ______          _____  
##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
##        Rsubread 1.16.1
## 
## //========================== subread-align setting ===========================\\
## ||                                                                            ||
## ||           Function : Read alignment                                        ||
## ||            Threads : 1                                                     ||
## ||         Input file : B_2.txt.gz                                            ||
## ||        Output file : B_2.bam (BAM)                                         ||
## ||         Index name : chr1                                                  ||
## ||       Phred offset : 33                                                    ||
## ||                                                                            ||
## ||          Min votes : 3                                                     ||
## ||         Max indels : 5                                                     ||
## ||  # of Best mapping : 1                                                     ||
## ||     Unique mapping : yes                                                   ||
## ||   Hamming distance : yes                                                   ||
## ||     Quality scores : no                                                    ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## //====================== Running (05-Mar-2015 12:27:06) ======================\\
## ||                                                                            ||
## || Decompress B_2.txt.gz...                                                   ||
## || The input file contains base space reads.                                  ||
## || Load the 1-th index block...                                               ||
## || Map reads...                                                               ||
## ||    0% completed,   0 mins elapsed, total=880k reads, rate=10.7k/s          ||
## ||    6% completed,   0 mins elapsed, total=880k reads, rate=10.8k/s          ||
## ||   12% completed,   0 mins elapsed, total=880k reads, rate=10.9k/s          ||
## ||   18% completed,   0 mins elapsed, total=880k reads, rate=10.9k/s          ||
## ||   24% completed,   0 mins elapsed, total=880k reads, rate=10.9k/s          ||
## ||   30% completed,   0 mins elapsed, total=880k reads, rate=10.9k/s          ||
## ||   36% completed,   0 mins elapsed, total=880k reads, rate=10.9k/s          ||
## ||   42% completed,   0 mins elapsed, total=880k reads, rate=10.7k/s          ||
## ||   48% completed,   0 mins elapsed, total=880k reads, rate=10.7k/s          ||
## ||   54% completed,   0 mins elapsed, total=880k reads, rate=10.6k/s          ||
## || Detect indels...                                                           ||
## ||   60% completed,   0 mins elapsed, total=880k reads, rate=10.0k/s          ||
## ||   62% completed,   0 mins elapsed, total=880k reads, rate=10.1k/s          ||
## ||   63% completed,   0 mins elapsed, total=880k reads, rate=10.2k/s          ||
## ||   64% completed,   0 mins elapsed, total=880k reads, rate=10.3k/s          ||
## ||   65% completed,   0 mins elapsed, total=880k reads, rate=10.4k/s          ||
## ||   66% completed,   0 mins elapsed, total=880k reads, rate=10.5k/s          ||
## ||   67% completed,   0 mins elapsed, total=880k reads, rate=10.6k/s          ||
## ||   68% completed,   0 mins elapsed, total=880k reads, rate=10.7k/s          ||
## ||   69% completed,   0 mins elapsed, total=880k reads, rate=10.8k/s          ||
## || Realign reads...                                                           ||
## ||   71% completed,   0 mins elapsed, total=880k reads, rate=10.5k/s          ||
## ||   72% completed,   1 mins elapsed, total=880k reads, rate=10.3k/s          ||
## ||   73% completed,   1 mins elapsed, total=880k reads, rate=10.0k/s          ||
## ||   74% completed,   1 mins elapsed, total=880k reads, rate=9.8k/s           ||
## ||   75% completed,   1 mins elapsed, total=880k reads, rate=9.6k/s           ||
## ||   76% completed,   1 mins elapsed, total=880k reads, rate=9.5k/s           ||
## ||   77% completed,   1 mins elapsed, total=880k reads, rate=9.3k/s           ||
## ||   78% completed,   1 mins elapsed, total=880k reads, rate=9.1k/s           ||
## ||   79% completed,   1 mins elapsed, total=880k reads, rate=9.0k/s           ||
## || 880316 reads were processed. Save the mapping results for them...          ||
## ||   81% completed,   1 mins elapsed, total=880k reads, rate=8.9k/s           ||
## ||   84% completed,   1 mins elapsed, total=880k reads, rate=8.9k/s           ||
## ||   86% completed,   1 mins elapsed, total=880k reads, rate=8.9k/s           ||
## ||   88% completed,   1 mins elapsed, total=880k reads, rate=8.9k/s           ||
## ||   90% completed,   1 mins elapsed, total=880k reads, rate=9.0k/s           ||
## ||   92% completed,   1 mins elapsed, total=880k reads, rate=9.0k/s           ||
## ||   94% completed,   1 mins elapsed, total=880k reads, rate=9.0k/s           ||
## ||   96% completed,   1 mins elapsed, total=880k reads, rate=9.0k/s           ||
## ||   98% completed,   1 mins elapsed, total=880k reads, rate=9.1k/s           ||
## ||                                                                            ||
## ||                          Completed successfully.                           ||
## ||                                                                            ||
## \\============================================================================//
## 
## //================================= Summary ==================================\\
## ||                                                                            ||
## ||          Processed : 880316 reads                                          ||
## ||             Mapped : 818107 reads (92.9%)                                  ||
## ||             Indels : 4829                                                  ||
## ||                                                                            ||
## ||       Running time : 1.6 minutes                                           ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
```

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:28:43 PST"
```

## Summarize mapped reads

Count numbers of reads mapped to NCBI Refseq genes.

> Summarize mapped reads to RefSeq genes. This will just take a few seconds. 
Note that the featureCounts function has built-in annotation for Refseq genes. 
featureCounts returns an R 'List' object, which includes raw read count for 
each gene in each library and also annotation information for genes such as 
gene identifiers and gene lengths. 


```r
setwd(datadir)
fc <- featureCounts(files=targets$OutputFile,annot.inbuilt="hg19")
```

```
## NCBI RefSeq annotation for hg19 (build 37.2) is used.
## 
##         ==========     _____ _    _ ____  _____  ______          _____  
##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
##        Rsubread 1.16.1
## 
## //========================== featureCounts setting ===========================\\
## ||                                                                            ||
## ||             Input files : 4 BAM files                                      ||
## ||                           S A_1.bam                                        ||
## ||                           S A_2.bam                                        ||
## ||                           S B_1.bam                                        ||
## ||                           S B_2.bam                                        ||
## ||                                                                            ||
## ||             Output file : ./.Rsubread_featureCounts_pid18594               ||
## ||             Annotations : /usr/local/lib/R/site-library/Rsubread/annot ... ||
## ||                                                                            ||
## ||                 Threads : 1                                                ||
## ||                   Level : meta-feature level                               ||
## ||              Paired-end : no                                               ||
## ||         Strand specific : no                                               ||
## ||      Multimapping reads : not counted                                      ||
## || Multi-overlapping reads : not counted                                      ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
## 
## //================================= Running ==================================\\
## ||                                                                            ||
## || Load annotation file /usr/local/lib/R/site-library/Rsubread/annot/hg19 ... ||
## ||    Features : 225074                                                       ||
## ||    Meta-features : 25702                                                   ||
## ||    Chromosomes : 52                                                        ||
## ||                                                                            ||
## || Process BAM file A_1.bam...                                                ||
## ||    Single-end reads are included.                                          ||
## ||    Assign reads to features...                                             ||
## ||    Total reads : 657474                                                    ||
## ||    Successfully assigned reads : 369440 (56.2%)                            ||
## ||    Running time : 0.04 minutes                                             ||
## ||                                                                            ||
## || Process BAM file A_2.bam...                                                ||
## ||    Single-end reads are included.                                          ||
## ||    Assign reads to features...                                             ||
## ||    Total reads : 562490                                                    ||
## ||    Successfully assigned reads : 312054 (55.5%)                            ||
## ||    Running time : 0.03 minutes                                             ||
## ||                                                                            ||
## || Process BAM file B_1.bam...                                                ||
## ||    Single-end reads are included.                                          ||
## ||    Assign reads to features...                                             ||
## ||    Total reads : 902841                                                    ||
## ||    Successfully assigned reads : 448593 (49.7%)                            ||
## ||    Running time : 0.05 minutes                                             ||
## ||                                                                            ||
## || Process BAM file B_2.bam...                                                ||
## ||    Single-end reads are included.                                          ||
## ||    Assign reads to features...                                             ||
## ||    Total reads : 880316                                                    ||
## ||    Successfully assigned reads : 435938 (49.5%)                            ||
## ||    Running time : 0.05 minutes                                             ||
## ||                                                                            ||
## ||                         Read assignment finished.                          ||
## ||                                                                            ||
## \\===================== http://subread.sourceforge.net/ ======================//
```

```r
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])
# Return to projdir
setwd(projdir)
print(Sys.time())
```

```
## [1] "2015-03-05 12:28:56 PST"
```

## Generate RPKM values

Generate RPKM values if you need them.


```r
x_rpkm <- rpkm(x,x$genes$Length)
print(Sys.time())
```

```
## [1] "2015-03-05 12:28:56 PST"
```

## Filter out low-count genes

> Only keep in the analysis those genes which had >10 reads per million mapped 
reads in at least two libraries. 


```r
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]
print(Sys.time())
```

```
## [1] "2015-03-05 12:28:56 PST"
```

## Perform voom normalization

> The figure below shows the mean-variance relationship estimated by voom for 
the data. 


```r
y <- voom(x,design,plot=TRUE)
```

![](rsubread_test_files/figure-html/unnamed-chunk-14-1.png) 

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:28:56 PST"
```

## Cluster libraries

> The following multi-dimensional scaling plot shows that sample A libraries 
are clearly separated from sample B libraries. 


```r
plotMDS(y,xlim=c(-2.5,2.5))
```

![](rsubread_test_files/figure-html/unnamed-chunk-15-1.png) 

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:28:56 PST"
```

## Assess differential expression

Fit linear model and assess differential expression.

> Fit linear models to genes and assess differential expression using the 
eBayes moderated t statistic. Here we list top 10 differentially expressed 
genes between B vs A. 


```r
fit <- eBayes(lmFit(y,design))
topTable(fit,coef=2)
```

```
##              GeneID Length logFC AveExpr   t P.Value adj.P.Val  B
## 100131754 100131754   1019   1.6      16 113 3.5e-28   6.3e-25 54
## 2023           2023   1812  -2.7      13 -91 2.2e-26   1.9e-23 51
## 2752           2752   4950   2.4      13  82 1.5e-25   9.1e-23 49
## 22883         22883   5192   2.3      12  64 1.8e-23   7.9e-21 44
## 6135           6135    609  -2.2      12 -62 3.1e-23   9.5e-21 44
## 6202           6202    705  -2.4      12 -62 3.2e-23   9.5e-21 44
## 4904           4904   1546  -3.0      11 -60 5.5e-23   1.4e-20 43
## 23154         23154   3705   3.7      11  55 2.9e-22   6.6e-20 41
## 8682           8682   2469   2.6      12  49 2.2e-21   4.3e-19 39
## 6125           6125   1031  -2.0      12 -48 3.1e-21   5.6e-19 39
```

```r
print(Sys.time())
```

```
## [1] "2015-03-05 12:28:57 PST"
```
