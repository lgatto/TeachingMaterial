TeachingMaterial
================

This repository is an aggregator for various [R](http://www.r-project.org/) teaching material.
Most of the courses are taught at the University of Cambridge, UK, and some have been adapted and exported outside. 
Each material subdirectory has its own repository; `TeachingMaterial` aggregates a snapshot as a central entry point. 
Aggregation is done using `git-subtree` (see the [administration page](https://github.com/lgatto/TeachingMaterial/wiki/TM-Administration) for details). 
The local copies linking to external repositories are prefixed with an underscore. 

Unless otherwise stated, all material is licensed under a 
[Creative Commons Attribution-ShareAlike 3.0 License](http://creativecommons.org/licenses/by-sa/3.0/). 
This means you are free to copy, distribute and transmit the work, 
adapt it to your needs as long as you cite its origin and, 
if you do redistribute it, do so under the same license.

If you like this material and/or this initiative, 
do not hesitate to let us know by starring the repo, 
tweeting about it and sharing it with your colleagues. 

## Material

### spr
- Description: Scientific Programmign with `R`, MPhil in Computational Biology
- Author: [Stephen Eglen](http://www.damtp.cam.ac.uk/user/sje30/) with contributions from [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/spr
- More details: https://github.com/lgatto/spr#readme

### QuickPackage
- Description: A very brief overview of `R` package creation
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/QuickPackage
- More details: https://github.com/lgatto/QuickPackage#readme

### R package development
- Description: Developing, documenting and testing an `R` package
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/RPackageDevelopment
- More details: https://github.com/lgatto/RPackageDevelopment#readme

### RBasics
- Description: A introduction to `R` for knowledgeable bioinformaticians. Used as `R` intro lecture for the [CSAMA](http://marray.economia.unimi.it/) workshop.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/RBasics
- More details: https://github.com/lgatto/RBasics#readme

### RIntro
- Description: A 1 day introductionary course, developed specifically for Microarray Analysis using `R` and Bioconductor.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/RIntro
- More details: https://github.com/lgatto/RIntro#readme

### R functional programming
- Description: A short topic section on functional programming in `R`
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/R-functional-programming
- More details: https://github.com/lgatto/R-functional-programming#readme

### R vectorisation
- Description: A short topic section on vectorisation in `R`
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: https://github.com/lgatto/R-vectorisation
- More details: https://github.com/lgatto/R-vectorisation#readme

### R debugging
- Description: A short topic section on debugging in `R`
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/R-debugging
- More details: https://github.com/lgatto/R-debugging#readme

### R parallel
- Description: A short topic section on parallal computing in `R`
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/R-parallel
- More details: https://github.com/lgatto/R-parallel#readme

### R object oriented programming
- Description: Covers S3, S4 and S4 Reference Classes OO programming using DNA/RNA sequence data manipulation as a working example.
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: https://github.com/lgatto/roo
- More details: https://github.com/lgatto/roo#readme

### Short S4 tutorial
- Description: A practical S4 object-oriented programming tutorial using microarrays as a working example.
- Author: [Laurent Gatto](https://github.com/lgatto) 
- Original repository: https://github.com/lgatto/S4-tutorial
- More details: https://github.com/lgatto/S4-tutorial#readme

### R and C/C++
- Description: Writing and calling `C`/`C++` code for/from `R`
- Author:  [Laurent Gatto](https://github.com/lgatto) 
- Many examples from the [Rcpp gallery](http://gallery.rcpp.org/) and 
  from [Dirk's slides](http://dirk.eddelbuettel.com/bio/presentations.html).
- Original repository: https://github.com/lgatto/rccpp
- More details: https://github.com/lgatto/rccpp#readme

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

