# R debugging and robust programming

25-26 February 2016 (Thursday-Friday), EMBL Heidelberg  
Instructors: Laurent Gatto, Robert Stojnic (University of Cambridge)  
Organiser: Wolfgang Huber (EMBL)  

This two-day course will teach participants debugging techniques and
good practice in writing reliable, robust code. The material will
provide the opportunity to gain experience and understanding of how to
identify, resolve, and avoid bugs, in order to produce
publication-quality code. The course will be taught using R and will
be driven by many practical exercises.  Course outline

The material will focus on:

- debugging, to fix problems with code
- defensive programming - writing effective tests to detect bugs
- profiling and optimisation of code

## Pre-requisites

The course is aimed at those with experience of scripting, who want to
learn more about writing robust and efficient code and who may want to
develop and release packages in the future.

## Content

#### [Part I:](https://github.com/lgatto/2016-02-25-adv-programming-EMBL/blob/master/01-intro.md)
- Coding style(s)
- Interactive use and programming
- Environments
- Tidy data
- Computing on the language

#### [Part II: Functional programming](https://github.com/lgatto/2016-02-25-adv-programming-EMBL/blob/master/02-funprog.md)

#### [Part III: Debugging](https://github.com/lgatto/2016-02-25-adv-programming-EMBL/blob/master/03-debug.md)
- Defensive programming
- Debbugging: techniques and tools
- Condition handling: try/tryCatch
- [Unit testing](https://github.com/lgatto/2016-02-25-adv-programming-EMBL/blob/master/unittesting.md)

#### [Part IV: Performance](https://github.com/lgatto/2016-02-25-adv-programming-EMBL/blob/master/04-perf.md)
- Benchmarking
- Profiling
- Optimisation
- Memory
- [Rcpp](https://github.com/lgatto/rccpp/blob/master/rc.md)

#### Other topics

- Packages and documentation
- Reproducible research and vignettes (`Rmarkdown`)
- Source code versioning with (for example) git and GitHub
- Automation with `Make`

## References

- [Previous courses](https://github.com/lgatto/teachingmaterial) and [here](https://github.com/DataProgrammers/2015-01-15-EMBLHeidelberg).
- [Advanced R](http://adv-r.had.co.nz/), Hadley Wickham.
- [The R Inferno](http://www.burns-stat.com/documents/books/the-r-inferno/), Patrick Burns.
- [An Introduction to the Interactive Debugging Tools in R](http://www.biostat.jhsph.edu/~rpeng/docs/R-debug-tools.pdf), Roger D. Peng.
- [R Programming for Bioinformatics](http://master.bioconductor.org/help/publications/books/r-programming-for-bioinformatics/), Robert Gentleman.


## License

This work is licensed under a CC BY-SA 3.0 License.
