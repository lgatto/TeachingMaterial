BIOSTAT 578A: Bioinformatics for Big Omics Data
===============================================

**Important Note:** I am in the process of modifying the content of this repository in preparation for Winter 2015. Stay tuned (I have now updated the first few lectures). If you want to be informed of all changes, please create a github account and watch the repository. Please also make sure you look at the "Getting Started" section below, as I expect you to do some things before the course actually starts. Please also login to `myuw` and look for other information on `canvas`.

**Instructor:** Raphael Gottardo, PhD, Fred Hutchinson Cancer Research Center

If you need to contact me, please email me at <rgottard@fhcrc.org>.

**Time and location:**
T & Th	9:00-10:20	HST T439

**Prerequisite:** BIOSTAT 511/12 or permission of the instructor. Please email me if you're unsure.

**Getting Started:** Please look at [this document](https://github.com/raphg/Biostat-578/blob/master/getting_started.md) to get you all set-up before the first class. This will include doing some reading/learning about R/Bioconductor/git/GitHub. 

**Grading scheme (Tentative):** HW (40%), Midterm (30%), Final project (30%)

**Important dates:** Midterm (Feb 19), Final project presentations: last 2 weeks of class (March 3, 5, 10 & 12).

**Scope:** This practical "hands-on" course in Bioinformatics for high dimensional omics will emphasize on how to use statistical methods, as well as the R programming language and the Bioconductor project, as tools to manipulate, visualize and analyze real world omics datasets. The course will be organized around the following topics:
- Introduction to computing for Bioinformatics using R: Introduction to R/RStudio, review of main data structures and tools for efficient and reproducible research, data manipulation and visualization
- Managing "big omics data" using relational databases: Overview of main database management systems (MySQL, Postgres, SQLite), and review of the Structured Query Language and main operations
- How to connect to a database from R, and alternative to databases in R (sqldf and data.table)
- How to evaluate and adjust the data for presence of "batch effect"
- Regression techniques for high throughput biomedical data: Multiple regression analysis and logistic regression, ANOVA and design of experiments
- Statistical methods for high dimensional hypothesis testing: Permutation tests, empirical Bayes and multiple comparison adjustment
- Modeling of gene expression data: Introduction to Bioconductor, and basic packages for gene expression analysis (GEOquery, Limma, DAVIDquery, etc)
- Genome-wide association studies and eQTLs; review of main packages in R/Bioconductor (e.g. rqtl)
- Overview of other high-throughput technologies (e.g. RNA-seq, ChIP-seq) and available tools in R/Bioconductor
- Data integration: Using R to integrate multiple data types and perform "systems biology" type analysis
- Drawbacks and limitations of high dimensional omics analysis (overfitting, inference)

*Note that this is tentative ouline and minor modifications are likely to occur. Please watch this page regularly for updates.*

**Lecture notes:**
Notes are provided as the source file (.Rmd) and resulting html file for online viewing. If you'd like to print these notes, please use the intermediate md file and [gitprint](http://gitprint.com/).
- 01/06/15 [Introduction to R](https://github.com/raphg/Biostat-578/blob/master/Introduction_to_R.Rmd) 
- 01/08/15 [Advanced graphics in R](https://github.com/raphg/Biostat-578/blob/master/Advanced_graphics_in_R.Rmd)
- 01/13/15 [Advanced data manipulation in R](https://github.com/raphg/Biostat-578/blob/master/Advanced_data_manipulation.Rmd)
- 01/15/15 Lab with Brian (in class). Use this opportunity to discuss HW1.
- 01/20/15 [Advanced data manipulation in R (suite)](https://github.com/raphg/Biostat-578/blob/master/Advanced_data_manipulation.Rmd) & [Molecular Biology 101](https://github.com/raphg/Biostat-578/blob/master/Biology_basics.Rmd)
