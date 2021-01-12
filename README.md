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

### Mass spectrometry and proteomics using R/Bioconductor
- Description: In this course, we will use R/Bioconductor packages to
  explore, process, visualise and understand mass spectrometry-based
  proteomics data, starting with raw data, and proceeding with
  identification and quantitation data, discussing some of their
  peculiarities compared to sequencing data along the way. The
  workflow is aimed at a beginner to intermediate level, such as, for
  example, seasoned R users who want to get started with mass
  spectrometry and proteomics, or proteomics practitioners who want to
  familiarise themselves with R and Bioconductor infrastructure.
- Direct link: http://bit.ly/bioc-ms-prot (see also this [3-days workshop](https://lgatto.github.io/2020-02-17-RProt-Berlin/))
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/bioc-ms-prot](https://github.com/lgatto/bioc-ms-prot)
- More details: [README](https://github.com/lgatto/bioc-ms-prot/blob/master/README.md)

### Visualising biomolecular data

- Description: This Visualisation of biomolecular data course is aimed
  at people who are already familiar with the R language and syntax,
  and who would like to get a hands-on introduction to visualisation,
  with a focus on biomolecular data in general, and proteomics in
  particular. This course is meant to be mostly hands-on, with an
  intuitive understanding of the underlying techniques.
- Direct link: http://bit.ly/biomolvis
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository:
  [https://github.com/lgatto/VisualisingBiomolecularData](https://github.com/lgatto/VisualisingBiomolecularData)

### A gentle introduction to git and Github
- Learning GitHub and a bit about git using a pancake recipe as an example.
- Blog post: https://lgatto.github.io/github-intro/
- Author: [Laurent Gatto](https://github.com/lgatto)

### Introduction to bioinformatics and data science
- Description: The WSBIM1207 course is an introduction to
  bioinformatics (and data science) for biology and biomedical
  students. It introduces bioinformatics methodology and technologies
  without relying on any prerequisites. The aim of this course is for
  students to be in a position to understand important notions of
  bioinformatics and tackle simple bioinformatics-related problems in
  R, in particular to develope simple R analysis scripts and
  reproducible analysis reports to interogate, visualise and
  understand data in a tidy tabular format.
- Direct link: [http://bit.ly/WSBIM1207](http://bit.ly/WSBIM1207)
- Author: [Laurent Gatto](https://github.com/lgatto)


### Bioinformatics
- Description: The
  [WSBIM1322](https://uclouvain.be/cours-2019-wsbim1322.html) course
  is teaches the basics of statistical data analysis applied to high
  throughput biology. It is aimed at biology and biomedical students
  that are already familiar with the R langauge (see the pre-requisits
  section below). The students will familiarise themselves with
  statitical learning concepts such as unsupervised and supervised
  learning, hypothesis testing, and extend their understanding and
  practive in R data structures and programming and the Bioconductor
  project.
- Direct link: [http://bit.ly/WSBIM1322](http://bit.ly/WSBIM1322)
- Author: [Laurent Gatto](https://github.com/lgatto)

### Advanced R programming

- Description: A two-day course taught on the 3-4 April 2017, teaching
  advanced techniques in writing reliable, robust code in R.
- Author: [Laurent Gatto](https://github.com/lgatto), and Robert
  Stojnic.
- Original repository:
  [https://github.com/lgatto/2017-04-03-adv-r-progr-EMBL](https://github.com/lgatto/2017-04-03-adv-r-progr-EMBL)
- Content: The material provides the opportunity to gain experience
  and understanding of object-oriented programming, packaging your
  code for distribution, advanced approaches for data visualisation,
  unit testing, and debugging.


### R debugging and robust programming

- Description: A 2-day workshop taught on the 25-26 February 2016 at
  the EMBL, Heidelberg. The course aims at teaching participants
  debugging techniques and good practice in writing reliable, robust
  code.
- Author: [Laurent Gatto](https://github.com/lgatto), based on
  previous content by Laurent Gatto and Robert Stojnic, and
  [*Advanced R*](http://adv-r.had.co.nz/), by Hadley Wickham.
- Original repository:
  [https://github.com/lgatto/2016-02-25-adv-programming-EMBL](https://github.com/lgatto/2016-02-25-adv-programming-EMBL)
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
- Original repository: [https://github.com/lgatto/rbc/](https://github.com/lgatto/rbc/)
- Content: `R` programming, plotting, `git`/`github` (via
  [software carpentry](http://software-carpentry.org/)), `make`,
  `shell` and `knitr`, profiling, testing, debugging.

### spr
- Description: Scientific Programming with `R`, MPhil in Computational Biology
- Author: [Stephen Eglen](http://www.damtp.cam.ac.uk/user/sje30/) with contributions from [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/spr](https://github.com/lgatto/spr)
- More details: [README](https://github.com/lgatto/spr#readme)

### Biostat-578
- Description: Bioinformatics for Big Omics Data
- Author: [Raphael Gottardo](http://www.rglab.org/), Fred Hutchinson Cancer Research Center
- Original repository: [https://github.com/raphg/Biostat-578](https://github.com/raphg/Biostat-578)
- More details: [README](https://github.com/raphg/Biostat-578/blob/master/README.md)

### github_tutorial
- Description: `git`/`github` guide
- Original repository: [https://github.com/kbroman/github_tutorial](https://github.com/kbroman/github_tutorial)
- Author: [Karl Broman](http://biostat.wisc.edu/~kbroman/), University of Wisconsin-Madison
- View it here: [http://kbroman.github.io/github_tutorial/](http://kbroman.github.io/github_tutorial/)

### minimal_make
- Description: minimal tutorial on `make`
- Original repository: [https://github.com/kbroman/minimal_make](https://github.com/kbroman/minimal_make)
- Author: [Karl Broman](http://biostat.wisc.edu/~kbroman/), University of Wisconsin-Madison
- View it here: [http://kbroman.github.io/minimal_make/](http://kbroman.github.io/minimal_make/)

### QuickPackage
- Description: Two brief overviews of `R` package creation
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/QuickPackage](https://github.com/lgatto/QuickPackage)
- More details: [README](https://github.com/lgatto/QuickPackage#readme)
- Download the [QuickPackage](https://github.com/lgatto/QuickPackage/blob/master/QuickPackage.pdf?raw=true) and [QuickPackageAndMore](https://github.com/lgatto/QuickPackage/blob/master/QuickPackageAndMore.pdf?raw=true) slides (pdf)

### R package development
- Description: Developing, documenting and testing an `R` package
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: [https://github.com/lgatto/RPackageDevelopment](https://github.com/lgatto/RPackageDevelopment)
- More details: [README](https://github.com/lgatto/RPackageDevelopment#readme)
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
- Original repository: [https://github.com/lgatto/RBasics](https://github.com/lgatto/RBasics)
- More details: [README](https://github.com/lgatto/RBasics#readme)
- Download the [pdf](https://github.com/lgatto/RBasics/blob/master/R-Basics.pdf?raw=true)

### RIntro
- Description: A 1 day introductionary course, developed specifically for Microarray Analysis using `R` and Bioconductor.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/RIntro](https://github.com/lgatto/RIntro)
- More details: [README](https://github.com/lgatto/RIntro#readme)
- Download the [pdf](https://github.com/lgatto/RIntro/blob/master/RIntro.pdf?raw=true)

### basicr
- Description: Material for the 2-day `R` introduction taught at
  Cambridge's
  [Graduate School of Life Sciences Training](http://www.training.cam.ac.uk/gsls/course/gsls-rintro).
- Authors: Robert Stojnic, [Laurent Gatto](https://github.com/lgatto),
  Rob Foy, David Molnar and John Davey, based on original slides by
  Ian Roberts and Robert Stojnic.
- Original repository: [https://github.com/johnomics/basicr/](https://github.com/johnomics/basicr/)
- Course page: [http://logic.sysbiol.cam.ac.uk/teaching/Rcourse/](http://logic.sysbiol.cam.ac.uk/teaching/Rcourse/)
- Download the slides for
  [day 1](https://github.com/johnomics/basicr/raw/master/Basic_R_Day_1_slides.pdf)
  and
  [day 2](https://github.com/johnomics/basicr/raw/master/Basic_R_Day_2_slides.pdf)
  in pdf.

### R functional programming
- Description: A short topic section on functional programming in `R`
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/R-functional-programming](https://github.com/lgatto/R-functional-programming)
- More details: [README](https://github.com/lgatto/R-functional-programming#readme)
- Download the [pdf](https://github.com/lgatto/R-functional-programming/blob/master/functional-programming.pdf?raw=true)

### R vectorisation
- Description: A short topic section on vectorisation in `R`
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/R-vectorisation](https://github.com/lgatto/R-vectorisation)
- More details: [README](https://github.com/lgatto/R-vectorisation#readme)
- Download the [pdf](https://github.com/lgatto/R-vectorisation/blob/master/vectorisation.pdf?raw=true)

### R debugging
- Description: A short topic section on debugging in `R`
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: [https://github.com/lgatto/R-debugging](https://github.com/lgatto/R-debugging)
- More details: [README](https://github.com/lgatto/R-debugging#readme)
- Download the [pdf](https://github.com/lgatto/R-debugging/blob/master/debugging.pdf?raw=true)

### R parallel
- Description: A short topic section on parallal computing in `R`
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: [https://github.com/lgatto/R-parallel](https://github.com/lgatto/R-parallel)
- More details: [README](https://github.com/lgatto/R-parallel#readme)
- Download the [pdf](https://github.com/lgatto/R-parallel/blob/master/parallel.pdf?raw=true)

### R object oriented programming
- Description: Covers S3, S4 and S4 Reference Classes OO programming using DNA/RNA sequence data manipulation as a working example.
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: [https://github.com/lgatto/roo](https://github.com/lgatto/roo)
- More details: [README](https://github.com/lgatto/roo#readme)
- Download the [pdf](https://github.com/lgatto/roo/blob/master/roo.pdf?raw=true)

### One day course on R OO programming and package development

- Description: A short 1-day course about R object-oriented
  programming, package development and various other topics (C
  interface, unit testing, debugging). See also other more recent and
  detailed lessons about these topic on this page.
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: [https://github.com/lgatto/advr1](https://github.com/lgatto/advr1)
- More details: [README](https://github.com/lgatto/advr1#readme)
- Download the [pdf](https://github.com/lgatto/advr1/blob/master/advr.pdf?raw=true)

### Short S4 tutorial
- Description: A practical S4 object-oriented programming tutorial using microarrays as a working example.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/S4-tutorial](https://github.com/lgatto/S4-tutorial)
- More details: [README](https://github.com/lgatto/S4-tutorial#readme)
- Download the [pdf](https://github.com/lgatto/S4-tutorial/blob/master/S4-tutorial.pdf?raw=true)

### R programming tutorial
- Description: A tutorial on R programming of intermediate level, focusing on some aspects of functional programming, profiling, testing, debugging and parallelisation. Used as more advanced `R` programming lecture during the [CSAMA](http://marray.economia.unimi.it/) workshop.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository: [https://github.com/lgatto/R-programming](https://github.com/lgatto/R-programming)
- More details: [README](https://github.com/lgatto/R-programming#readme)
- Download the [pdf](https://github.com/lgatto/R-programming/blob/master/R-programming.pdf?raw=true)

### R and C/C++
- Description: Writing and calling `C`/`C++` code for/from `R`
- Author:  [Laurent Gatto](https://github.com/lgatto)
- Many examples from the [Rcpp gallery](http://gallery.rcpp.org/) and
  from [Dirk's slides](http://dirk.eddelbuettel.com/bio/presentations.html).
- Original repository: [https://github.com/lgatto/rccpp](https://github.com/lgatto/rccpp)
- More details: [README](https://github.com/lgatto/rccpp#readme)
- Download the [pdf](https://github.com/lgatto/rccpp/blob/master/rccpp.pdf?raw=true)

### visualisation
- Description: visualising data using `ggplot2` (more to come soon)
- Author: Mark Dunning
- Original repository:
  [https://github.com/lgatto/visualisation](https://github.com/lgatto/visualisation)
- Download the [pdf](https://github.com/lgatto/visualisation/blob/master/ggplot2_cambr28oct2013.pdf?raw=true)

### sequences
- Description: Educational package used in `R` to illustrate OO programming and package development
- Author: [Laurent Gatto](https://github.com/lgatto) and Robert Stojnić
- Original repository: [https://github.com/lgatto/sequences](https://github.com/lgatto/sequences)
- More details: [DESCRIPTION](https://github.com/lgatto/sequences/blob/master/DESCRIPTION)
- Installation from CRAN: `install.packages("sequences")`
- Installation from github (requires `R` and `C/C++` building tools):

```c
library(devtools)
install_github("lgatto/sequences")
```

### Best practices in bioinformatics research: open source software and reproducibility

- Description: A short course for a Bioinformatics minor at the
  University of Cambridge. What is open science (data, source/code,
  access), and how can we enable it? What is reproducible research,
  and why do we need it and how can we implement it? The objective is
  to familiarise students with concepts and tools of open science and
  reproducible research.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository:
  [https://github.com/lgatto/open-rr-bioinfo-best-practice](https://github.com/lgatto/open-rr-bioinfo-best-practice)

### Beginner's statistics in R

- Description: A 2.5 days introductionary course focusing on R and
  basics statistics for proteomics scientists. The R intruduction
  material is based on the
  [Data Carpentry R analysis lesson](http://www.datacarpentry.org/R-ecology-lesson/)
  and leads to the introduction and application of basic uni-variate
  statistics using proteomics data. The course was developed and
  taught as part of the May Institute, at the Northeastern University,
  Boston, MA in May 2017.
- Authors: [Laurent Gatto](https://github.com/lgatto) and
  [Meena Choi](https://github.com/MeenaChoi), with material from the Data Carpentry R
  lesson.
- Original repository:
  [https://github.com/lgatto/2017-05-03-RstatsIntro-NEU](https://github.com/lgatto/2017-05-03-RstatsIntro-NEU)

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
- Original repository:
  [https://github.com/lgatto/Quantitative-Proteomics-and-Data-Analysis](https://github.com/lgatto/Quantitative-Proteomics-and-Data-Analysis)
- Slides (http://bit.ly/qprotda) and reproducible vignette (http://bit.ly/qprotdavig)

### R and Bioconductor for Mass Spectrometry and Proteomics data analysis

- Description: Material for the 2-day R/Bioconductor Proteomics
  Workshop at Stellenbosch University, October 2016
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository:
  [https://github.com/lgatto/2016-10-rbioc-prot-Stellenbosch](https://github.com/lgatto/2016-10-rbioc-prot-Stellenbosch)

### An Introduction to Machine Learning with R

- Description: This introductory workshop on machine learning with R
  is aimed at participants who are not experts in machine learning
  (introductory material will be presented as part of the course), but
  have some familiarity with scripting in general and R in
  particular. The workshop will offer a hands-on overview of typical
  machine learning applications in R, including unsupervised
  (clustering, such as hierarchical and k-means clustering, and
  dimensionality reduction, such as principal component analysis) and
  supervised (classification and regression, such as K-nearest
  neighbour and linear regression) methods. We will also address
  questions such as model selection using cross-validation.
- Author: [Laurent Gatto](https://github.com/lgatto)
- Original repository:
  [https://lgatto.github.io/IntroMachineLearningWithR/](https://lgatto.github.io/IntroMachineLearningWithR/)
- Direct access to the material: [bookdown formatted](https://lgatto.github.io/IntroMachineLearningWithR/)

## License

We try to only aggregate material that is openly available, generally
under
[Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/),
which gives you the right to share and adapt the material as long as
you credit to original author(s). Please refer to the orignal
repository for details.
