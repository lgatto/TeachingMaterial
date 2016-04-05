TeachingMaterial
================

This repository is an aggregator for various
[R](http://www.r-project.org/), `make` and `git`/`github` teaching
material.  Most of the courses are taught at the University of
Cambridge, UK, and some have been adapted and exported outside. We
would also like to acknowledge contributions from
[Aleksandra Pawlik](http://www.software.ac.uk/about/people/aleksandra-pawlik),
Software Sustainability Institute,
[Raphael Gottardo](http://www.rglab.org/), Fred Hutchinson Cancer
Research Center and [Karl Broman](http://biostat.wisc.edu/~kbroman/),
University of Wisconsin-Madison.

Each material subdirectory has its own repository; `TeachingMaterial`
aggregates a snapshot as a central entry point.  Aggregation is done
using `git-subtree` (see the
[administration page](https://github.com/lgatto/TeachingMaterial/wiki/TM-Administration)
for details).  The local copies linking to external repositories are
prefixed with an underscore.

Unless otherwise stated, all material is licensed under a
[Creative Commons Attribution-ShareAlike 3.0 License](http://creativecommons.org/licenses/by-sa/3.0/).
This means you are free to copy, distribute and transmit the work,
adapt it to your needs as long as you cite its origin and, if you do
redistribute it, do so under the same license.

See also the
[`TeachingMaterial` wiki](https://github.com/lgatto/teachingmaterial/wiki)
for meta-information about the repository and general `R` installation
material and links.

If you like this material and/or this initiative, do not hesitate to
let us know by starring the repo, tweeting about it and sharing it
with your colleagues.

## Material

### R debugging and robust programming

- Description: A 2-day workshop taught on the 25-26 February 2016 at
  the EMBL, Heidelberg. The course aims at teaching participants
  debugging techniques and good practice in writing reliable, robust
  code.
- Author: [Laurent Gatto](https://github.com/lgatto), based on
  previous content by Laurent Gatto and Robert Stojnic, and
  [*Advanced R*](http://adv-r.had.co.nz/), by Hadley Wickham.
- Original repository: https://github.com/lgatto/2016-02-25-adv-programming-EMBL
- Content: Part I: Coding style(s), Interactive use and programming,
  Environments, Tidy data, Computing on the language. Part II:
  Functions, Robust programming with functions, Scoping, Closures,
  High-level functions, Vectorisation. Part III: Defensive
  programming, Debbugging: techniques and tools, Condition handling:
  try/tryCatch, Unit testing. Part IV: Benchmarking, Profiling,
  Optimisation, Memory, Rcpp.
- More details: https://github.com/lgatto/2016-02-25-adv-programming-EMBL/blob/master/README.md

### rbc
- Description:
  [Software carpentry R bootcamp](http://github.com/lgatto/rbc), Jan
  7-8, 2014, Cambridge, UK and 6-7 Nov 2014, Zurich, Switzerland.
- Authors: [Stephen Eglen](http://www.damtp.cam.ac.uk/user/sje30/),
  [Laurent Gatto](https://github.com/lgatto), Robert Stojnić and
  [Aleksandra Pawlik](http://www.software.ac.uk/about/people/aleksandra-pawlik)
- Original repository: https://github.com/lgatto/rbc/
- Content: `R` programming, plotting, `git`/`github` (via
  [software carpentry](http://software-carpentry.org/)), `make`,
  `shell` and `knitr`, profiling, testing, debugging.

### spr
- Description: Scientific Programming with `R`, MPhil in Computational Biology
- Author: [Stephen Eglen](http://www.damtp.cam.ac.uk/user/sje30/) with contributions from [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/spr
- More details: https://github.com/lgatto/spr#readme

### Biostat-578
- Description: Bioinformatics for Big Omics Data
- Author: [Raphael Gottardo](http://www.rglab.org/), Fred Hutchinson Cancer Research Center
- Original repository: https://github.com/raphg/Biostat-578
- More details: https://github.com/raphg/Biostat-578/blob/master/README.md

### rbioc-proteomics
- Using R/Bioconductor for MS-based proteomics data analysis
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/rbioc-proteomics
- More details: https://github.com/lgatto/rbioc-proteomics#readme

### github_tutorial
- Description: `git`/`github` guide
- Original repository: https://github.com/kbroman/github_tutorial
- Author: [Karl Broman](http://biostat.wisc.edu/~kbroman/), University of Wisconsin-Madison
- View it here: http://kbroman.github.io/github_tutorial/

### minimal_make
- Description: minimal tutorial on `make`
- Original repository: https://github.com/kbroman/minimal_make
- Author: [Karl Broman](http://biostat.wisc.edu/~kbroman/), University of Wisconsin-Madison
- View it here: http://kbroman.github.io/minimal_make/

### QuickPackage
- Description: Two brief overviews of `R` package creation
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/QuickPackage
- More details: https://github.com/lgatto/QuickPackage#readme
- Download the [QuickPackage](https://github.com/lgatto/QuickPackage/blob/master/QuickPackage.pdf?raw=true) and [QuickPackageAndMore](https://github.com/lgatto/QuickPackage/blob/master/QuickPackageAndMore.pdf?raw=true) slides (pdf)

### R package development
- Description: Developing, documenting and testing an `R` package
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/RPackageDevelopment
- More details: https://github.com/lgatto/RPackageDevelopment#readme
- Download the [pdf](https://github.com/lgatto/RPackageDevelopment/blob/master/rpd.pdf?raw=true)

### Benchmarking, profiling and optimisation
- Description: Benchmarking, profiling and optimisation
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/R-bmark-prof-optim
- More details: https://github.com/lgatto/R-bmark-prof-optim#readme
- Read the [material](https://github.com/lgatto/R-bmark-prof-optim/blob/master/bmark-prof-optim.md)

### RBasics
- Description: A introduction to `R` for knowledgeable bioinformaticians. Used as `R` intro lecture for the [CSAMA](http://marray.economia.unimi.it/) workshop.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/RBasics
- More details: https://github.com/lgatto/RBasics#readme
- Download the [pdf](https://github.com/lgatto/RBasics/blob/master/R-Basics.pdf?raw=true)

### RIntro
- Description: A 1 day introductionary course, developed specifically for Microarray Analysis using `R` and Bioconductor.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/RIntro
- More details: https://github.com/lgatto/RIntro#readme
- Download the [pdf](https://github.com/lgatto/RIntro/blob/master/RIntro.pdf?raw=true)

### basicr
- Description: Material for the 2-day `R` introduction taught at Cambridge's [Graduate School of Life Sciences Training](http://www.training.cam.ac.uk/gsls/course/gsls-rintro).
- Authors: Robert Stojnic, [Laurent Gatto](https://github.com/lgatto), Rob Foy, David Molnar and John Davey, based on original slides by Ian Roberts and Robert Stojnic.
- Original repository: https://github.com/johnomics/basicr/
- Course page: http://logic.sysbiol.cam.ac.uk/teaching/Rcourse/
- Download the slides for [day 1](https://github.com/johnomics/basicr/raw/master/Basic_R_Day_1_slides.pdf) and [day 2](https://github.com/johnomics/basicr/raw/master/Basic_R_Day_2_slides.pdf) in pdf.

### R functional programming
- Description: A short topic section on functional programming in `R`
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/R-functional-programming
- More details: https://github.com/lgatto/R-functional-programming#readme
- Download the [pdf](https://github.com/lgatto/R-functional-programming/blob/master/functional-programming.pdf?raw=true)

### R vectorisation
- Description: A short topic section on vectorisation in `R`
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/R-vectorisation
- More details: https://github.com/lgatto/R-vectorisation#readme
- Download the [pdf](https://github.com/lgatto/R-vectorisation/blob/master/vectorisation.pdf?raw=true)

### R debugging
- Description: A short topic section on debugging in `R`
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/R-debugging
- More details: https://github.com/lgatto/R-debugging#readme
- Download the [pdf](https://github.com/lgatto/R-debugging/blob/master/debugging.pdf?raw=true)

### R parallel
- Description: A short topic section on parallal computing in `R`
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/R-parallel
- More details: https://github.com/lgatto/R-parallel#readme
- Download the [pdf](https://github.com/lgatto/R-parallel/blob/master/parallel.pdf?raw=true)

### R object oriented programming
- Description: Covers S3, S4 and S4 Reference Classes OO programming using DNA/RNA sequence data manipulation as a working example.
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/roo
- More details: https://github.com/lgatto/roo#readme
- Download the [pdf](https://github.com/lgatto/roo/blob/master/roo.pdf?raw=true)

### Short S4 tutorial
- Description: A practical S4 object-oriented programming tutorial using microarrays as a working example.
- Author: [Laurent Gatto](https://github.com/lgatto) 
- Original repository: https://github.com/lgatto/S4-tutorial
- More details: https://github.com/lgatto/S4-tutorial#readme
- Download the [pdf](https://github.com/lgatto/S4-tutorial/blob/master/S4-tutorial.pdf?raw=true)

### R programming tutorial
- Description: A tutorial on R programming of intermediate level, focusing on some aspects of functional programming, profiling, testing, debugging and parallelisation. Used as more advanced `R` programming lecture during the [CSAMA](http://marray.economia.unimi.it/) workshop.
- Author: [Laurent Gatto](https://github.com/lgatto) 
- Original repository: https://github.com/lgatto/R-programming
- More details: https://github.com/lgatto/R-programming#readme
- Download the [pdf](https://github.com/lgatto/R-programming/blob/master/R-programming.pdf?raw=true)

### R and C/C++
- Description: Writing and calling `C`/`C++` code for/from `R`
- Author:  [Laurent Gatto](https://github.com/lgatto) 
- Many examples from the [Rcpp gallery](http://gallery.rcpp.org/) and 
  from [Dirk's slides](http://dirk.eddelbuettel.com/bio/presentations.html).
- Original repository: https://github.com/lgatto/rccpp
- More details: https://github.com/lgatto/rccpp#readme
- Download the [pdf](https://github.com/lgatto/rccpp/blob/master/rccpp.pdf?raw=true)

### visualisation
- Description: visualising data using `ggplot2` (more to come soon)
- Author: Mark Dunning
- Original repository: https://github.com/lgatto/visualisation
- Download the [pdf](https://github.com/lgatto/visualisation/blob/master/ggplot2_cambr28oct2013.pdf?raw=true)

### sequences
- Description: Educational package used in `R` to illustrate OO programming and package development
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/sequences
- More details: https://github.com/lgatto/sequences/blob/master/DESCRIPTION
- Installation from CRAN: `install.packages("sequences")`
- Installation from github (requires `R` and `C/C++` building tools): 

```c
library(devtools)
install_github("sequences", "lgatto")
```

### Best practices in bioinformatics research: open source software and reproducibility

- Description: A short course for a Bioinformatics minor at the
  University of Cambridge. What is open science (data, source/code,
  access), and how can we enable it? What is reproducible research,
  and why do we need it and how can we implement it? The objective is
  to familiarise students with concepts and tools of open science and
  reproducible research.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/open-rr-bioinfo-best-practice

### Statistics primer

- Description: A short course for a Bioinformatics minor at the
  University of Cambridge. Introducing basic concepts in statistics:
  experimental design, randomisation, technical and biological
  variation, power analysis, hypothesis testing, confidence interval,
  what is a p-value, false discovery rate, multiple testing
  adjustment, dangers of uninformed statistical practice. The objective
  is familiarise students with basic statistical concepts and initiate
  them to statistical thinking. They should be able to critically
  assess an experimental design and the reporting of a simple
  statistical analysis.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/statistics-primer
- Slides on
  [experimental desing](https://github.com/lgatto/statistics-primer/blob/master/expdes-slides.pdf),
  [significance testing](https://github.com/lgatto/statistics-primer/blob/master/test-slides.pdf), 
  and
  [practical](https://htmlpreview.github.io/?https://github.com/lgatto/statistics-primer/blob/master/03-practical.html).

### Inspection, visualisation and analysis of quantitative proteomics data

- Description: A discussion I lead in the frame of the
[Quantitative Proteomics and Data Analysis](https://www.biochemistry.org/Events/tabid/379/View/Programme/MeetingNo/TD007/Default.aspx)
Course.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/Quantitative-Proteomics-and-Data-Analysis
- Slides (http://bit.ly/qprotda) and reproducible vignette (http://bit.ly/qprotdavig)
