Bioinformatics for Big Omics Data: Introduction to R
========================================================
width: 1440
height: 900
transition: none
font-family: 'Helvetica'
css: my_style.css
author: Raphael Gottardo, PhD
date: January 09, 2014

<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png" /></a><br /><tiny>This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US">Creative Commons Attribution-ShareAlike 3.0 Unported License</tiny></a>.

What is Bioinformatics?
=======================

- It depends who you ask but, according to [Wikipedia](http://www.wikipedia.org): 


> Bioinformatics is an interdisciplinary field that develops and improves on methods for storing, retrieving, organizing and analyzing biological data.

- Also according to [Wikipedia](http://www.wikipedia.org): 

> Bioinformatics uses many areas of computer science, mathematics and engineering to process biological data. 

- **What about statistics?** Statistics plays an important role in Bioinformatics. You can't analyze data without statistics.

- **Bioinformatics is an art** that involves biology, statistics, mathematics, programming, and discipline (e.g. reproducibility).


Outline
=======

The syllabus is available on [GitHub](https://github.com/) at [https://github.com/raphg/Biostat-578](https://github.com/raphg/Biostat-578)

- All lecture notes, R code, etc, are (will) also be available there.

    - Lecture notes are written in R markdown and processed using knitr, and thus are **fully reproducible** including R code and examples. The notes are fully versioned via git (and [GitHub](https://github.com/)). 

    - You will learn more about these tools later. 
    
- You are more than welcome to modify, correct and/or contribute. Very easy to do via **pull requests**. Please do not send me correction via emails. 

- **Grading scheme (Tentative)** HW: 40%, Midterm: 30%, Final project: 30%


Reproducible research
=====================

Throughout this course, I will emphasize the use of open source tools for reproducible  research. You might want to read the following article:

Huang, Y., & Gottardo, R. (2012). Comparability and reproducibility of biomedical data. *Briefings in bioinformatics*. doi:10.1093/bib/bbs078

I expect each one of you to become familiar with:

- [R](r-project.org) and [RStudio](r-project.org)

- [knitr](http://yihui.name/knitr/). An R package to provide dynamic report generation with R. 

- [git](http://git-scm.com/) and [GitHub](https://github.com/). RStudio provides an integrated client for git. 

We will make heavy use of these tools throughout the course. 

Some R history
========================================================

- R is the son of S
- S is a statistical programming language developed by John Chambers from Bell Labs
- Goal of S was "to turn ideas into software, quickly and faithfully"
- S was created in 1976
- New S language arrived in 1988 (Blue Book) and introduced many changes (macros to functions)
- Version 4 was introduced in 1998 and introduced a formal class-method model
- In 1993, StatSci (maker of S-Plus) acquire exclusive license to S
- S-Plus integrates S with a nice GUI interface and full customer support
- R was created by Ross Ihaka and Robert Gentleman at the University of Auckland, New Zealand

R
===
- The R project started in 1991
- R first appeared in 1996 as an open-source software!
- **Highly customizable via packages**
- R based community, power of collaboration with thousands of packages freely available
- Many of my favorite R capabilities are not part of the base distribution
- Many graphical user interface to R both free and commercial (e.g. R studio and Revolution)

What is R?
==========
R is an integrated suite of software facilities for data manipulation, calculation and graphical display. It includes:
- an effective data handling and storage facility
- a suite of operators for calculations on arrays, in particular matrices
- a large, coherent, integrated collection of intermediate tools for data analysis
- graphical facilities for data analysis and display either on-screen or on hardcopy, and
- well-developed, simple and effective programming language which includes conditionals, loops, user-defined recursive functions and input and output facilities.

R in the NY Times
=================

"Despite" being free and open-source, R is widely used by data analysts inside corporations and academia.

See [NY Times](http://www.nytimes.com/2009/01/07/technology/business-computing/07program.html?pagewanted=all&_r=0) article

![R in the NY Times](Introduction_to_R-figure/R_NY_Times.png)


Some references
===============

Some references to get you started if you need to brush up your R skills.
- [aRrgh: a newcomer's (angry) guide to R](http://tim-smith.us/arrgh/) by Tim Smith and Kevin Ushey
- [Introductory Statistics with ](http://www.amazon.com/Introductory-Statistics-R-Computing/dp/0387790535) by Peter Dalgaard
- R reference card http://cran.r-project.org/doc/contrib/Short-refcard.pdf
- R tutorial http://www.cyclismo.org/tutorial/R/
- [R project](http://r-project.org) and [Bioconductor](bioconductor.org)

More advanced:

- [Hadley Wickham's book](http://adv-r.had.co.nz/)


RStudio
=======

[RStudio](http://www.rstudio.com/) is a free and open source integrated development environment. 

![RStudio IDE](http://www.rstudio.com/images/screenshots/rstudio-windows.png)

-------------
- Cross platform
- Syntax highlighting, code completion, and smart indentation
- Execute R code directly from the source editor
- Easily manage multiple working directories using projects
- Plot history, zooming, and flexible image and PDF export
- Integrated with knitr
- Integrated with Git for version control


R basics
========
R is an overgrown calculator!

```r
2+2
```

```
[1] 4
```

```r
exp(-2)
```

```
[1] 0.1353
```

```r
pi
```

```
[1] 3.142
```

```r
sin(2*pi)
```

```
[1] -2.449e-16
```


Getting help
============
You can easily find help via the command line:


```r
help(pi) # equivalent ?pi
?sqrt
?sin
?Special
```


If you don't know the exact name, use

```r
help.search("trigonometry")
??trigonometry
```

Or using the help tab integrated in RStudio, or using your favorite search engine!

Assignment
==========

Need a way to store intermediate results:

```r
x <- 2
y <- 2
x+y
```

```
[1] 4
```


Try to use meaningful names! 
Have a look at:
- [Hadley Wickham's book](http://adv-r.had.co.nz/)
- [Google's coding standards](http://google-styleguide.googlecode.com/svn/trunk/Rguide.xml)


Vectorized arithmetic
========================================================
We cannot do much statistics with a single number!

We need a way to store a sequence/list of numbers

One can simply concatenate elements with the `c` function.


```r
weight <- c(60, 72, 75, 90, 95, 72)
weight[1]
```

```
[1] 60
```

```r
height <- c(1.75, 1.80, 1.65, 1.90, 1.74, 1.91)
bmi <- weight/height^2 # vector based operation
bmi
```

```
[1] 19.59 22.22 27.55 24.93 31.38 19.74
```

- Vector based operation are much faster!
- `c` can be used to concatenate strings and numbers.

**Exercise:** Find at least one other way to create a vector.

Data structures
========================================================

Even vectors can be limited and we need richer structures.

**Homogeneous:**
- Vectors (1-d)
- Matrix (2-d)
- Arrays (n-d)

Can be logical, integer, double (often called numeric), or character

**Heterogeneous:**
- List
- Dataframes

For more details: [http://adv-r.had.co.nz/Data-structures.html](http://adv-r.had.co.nz/Data-structures.html)

Vectors
=======

We have three types of vectors: numeric, logical, character

```r
# Numeric vectors
x <- c(1, 5, 8)
x
```

```
[1] 1 5 8
```

```r
# Logical vectors
x <- c(TRUE, TRUE, FALSE, TRUE)
x
```

```
[1]  TRUE  TRUE FALSE  TRUE
```

```r
# Character vectors
x <- c("Hello", "my", "name", "is", "Francis")
x
```

```
[1] "Hello"   "my"      "name"    "is"      "Francis"
```


**Exercise:** Create a vector with the following elements 1, 3, 10, -1, call your vector x. Take the square root of x. Take the log of (1+x).

Missing and special values
==========================
We have already encountered the `NaN` symbol meaning not-a-number, and `Inf`, `-Inf`. In practical data analysis a data point is frequently unavailable. In R, missing values are denoted by `NA`. 

Depending on the context, R provides different ways to deal with missing values.


```r
weight <- c(60, 72, 75, 90, NA, 72)
mean(weight)
```

```
[1] NA
```

```r
mean(weight, na.rm=TRUE)
```

```
[1] 73.8
```



Matrices and arrays
===================

A matrix is a two dimensional array of numbers. Matrices can be used to perform statistical operations (linear algebra). However, they can also be used to hold tables. 


```r
x <- 1:12
length(x)
```

```
[1] 12
```

```r
dim(x)
```

```
NULL
```

```r
dim(x) <- c(3, 4)
x
```

```
     [,1] [,2] [,3] [,4]
[1,]    1    4    7   10
[2,]    2    5    8   11
[3,]    3    6    9   12
```

```r
x <- matrix(1:12, nrow=3, byrow=TRUE)
x <- matrix(1:12, nrow=3, byrow=FALSE)
rownames(x) <- c("A", "B", "C")
colnames(x) <- c("1", "2", "x", "y")
```


Matrices and Arrays (suite)
===================

Matrices can also be formed by "glueing" rows and columns using `cbind` and `rbind`. This is the equivalent of `c` for vectors.


```r
x1 <- 1:4
x2 <- 5:8
y1 <- c(3, 9)
my_matrix <- rbind(x1, x2)
new_matrix <- cbind(my_matrix, y1)
new_matrix
```

```
           y1
x1 1 2 3 4  3
x2 5 6 7 8  9
```


n-dimesional arrays generalize matrices, as follows:

```r
array(1:9, c(3, 3, 3))
```


Factors
=======

It is common to have categorical data in statistical data analysis (e.g. Male/Female). In R such variables are referred to as factors. Makes it possible to assign meaningful names to categories. A factor has a set of levels.


```r
pain <- c(0, 3, 2, 2, 1)
fpain <- factor(pain)
levels(fpain) <- c("none", "mild", "medium", "severe")
is.factor(fpain)
```

```
[1] TRUE
```

```r
is.vector(fpain)
```

```
[1] FALSE
```

```r
# Additional attribute
levels(fpain)
```

```
[1] "none"   "mild"   "medium" "severe"
```

A factor is very similar to an integer vector with a set of labels. While factors look like character vectors, they are not. So be careful when converting factors to characters and vice-versa. For example, use `stringsAsFactors = FALSE` when reading dataframes (more on this later).

Lists
=====
Lists can be used to store objects (of possibly different kinds/sizes) into a larger composite object. 

```r
x <- c(31, 32, 40)
y <- factor(c("F", "M", "M", "F"))
# Different types and dimensions!
z <- c("London", "School")
my_list <- list(age=x, sex=y, meta=z)
my_list
```

```
$age
[1] 31 32 40

$sex
[1] F M M F
Levels: F M

$meta
[1] "London" "School"
```

```r
my_list$age
```

```
[1] 31 32 40
```


Data Frames
==========
A data frame is a "data matrix" or a "data set". It is a list of vectors and/or factors of the same length that are related "across" such that data in the same position come from the same experimental unit (subject, gene, etc).

```r
my_df <- data.frame(age=c(31, 32, 40, 50), sex=c("M", "M", "F", "M"))
my_df$age
```

```
[1] 31 32 40 50
```

Why do we need data frames if it is simply a list? 
- More efficient storage, and indexing!

Dataframes are similar to database tables.
R provides some (more efficient) alternatives to dataframe. More later!

Names
==========
Name(s) of an R object can be accessed and/or modified with the `names` function (method).


```r
x <- rep(1:3)
names(x)
```

```
NULL
```

```r
names(x) <- c("a", "b", "c")
my_df <- data.frame(age=c(31,32,40,50), sex=y)
my_df
```

```
  age sex
1  31   F
2  32   M
3  40   M
4  50   F
```

```r
names(my_df)
```

```
[1] "age" "sex"
```

```r
names(my_df) <- c("age", "gender")
names(my_df)[1] <- c("Age")
```


Names are a special kind of attributes. See more here: http://adv-r.had.co.nz/Data-structures.html#attributes


Indexing
========

Indexing is a great way to directly 
access elements of interest.


```r
# Indexing a vector
pain <- c(0, 3, 2, 2, 1)
pain[1]
```

```
[1] 0
```

```r
pain[2]
```

```
[1] 3
```

```r
pain[1:2]
```

```
[1] 0 3
```

```r
pain[c(1, 3)]
```

```
[1] 0 2
```

```r
pain[-5]
```

```
[1] 0 3 2 2
```

Note that with a data frame, the indexing of subject is straightforward!

Indexing (suite)
================


```r
# Indexing a matrix
my_matrix[1, 1]
```

```
x1 
 1 
```

```r
my_matrix[1, ]
```

```
[1] 1 2 3 4
```

```r
my_matrix[, 1]
```

```
x1 x2 
 1  5 
```

```r
my_matrix[, -2]
```

```
   [,1] [,2] [,3]
x1    1    3    4
x2    5    7    8
```

```r
# Indexing list is done in the same way
my_list[3]
```

```
$meta
[1] "London" "School"
```

```r
my_list[[3]]
```

```
[1] "London" "School"
```

```r
my_list[[3]][1]
```

```
[1] "London"
```

```r
# Indexing a data frame
my_df[1, ]
```

```
  Age gender
1  31      F
```

```r
my_df[2, ]
```

```
  Age gender
2  32      M
```




Indexing by name
================


```r
my_list$age
```

```
[1] 31 32 40
```

```r
my_list["age"]
```

```
$age
[1] 31 32 40
```

```r
my_list[["age"]]
```

```
[1] 31 32 40
```

```r
my_df["Age"]
```

```
  Age
1  31
2  32
3  40
4  50
```

```r
# Try also
# my_df[1]
# my_df[[1]]
```


What is the main difference between `[[]]` and `[]`?

Conditional indexing
====================
Indexing can be conditional on another variable!


```r
pain <- c(0, 3, 2, 2, 1)
sex <- factor(c("M", "M", "F", "F", "M"))
age <- c(45, 51, 45, 32, 90)
pain[sex=="M"]
```

```
[1] 0 3 1
```

```r
pain[age>32]
```

```
[1] 0 3 2 1
```


**Exercise:** Do the same by indexing with F.
Do the same with age less than 80.

Functions and arguments
=======================

Many things in R are done using function calls, commands that look like an application of a mathematical function of one or several variables, e.g. `log(x)`, `plot(weight,height)`

We will see more on this when explore advance graphics in R.

Most function arguments have sensible default and can thus be omitted, e.g. `plot(weight, height, col=1)`

**Note:** If you do not specify the names of the argument, the order is important!

Loops and conditional statements
================================

R is a true programming language, and thus has a rich syntax including `for` loops and conditional statements (`while`, `if`, `ifelse`, etc).


```r
# A simple if statement
x <- -2
if(x>0) {
  print(x)
} else {
  print(-x)
}
```

```
[1] 2
```

```r

if(x>0) {
  print(x)
} else if(x==0) {
  print(0)
} else {
  print(-x)
}
```

```
[1] 2
```


--------------


```r
# For loops
n <- 1000000
x <- rnorm(n, 10, 1)
y <- x^2
y <- rep(0, n)
for(i in 1:n) {
  y[i]<-sqrt(x[i])
}

y[1:10]
```

```
 [1] 3.234 3.515 3.110 3.215 3.031 3.292 3.426 3.405 3.356 2.993
```

```r

# While loops
counter <- 1
while(counter<=n) {
  y[counter] <- sqrt(x[counter])
  counter <- counter+1
}

y[1:10]
```

```
 [1] 3.234 3.515 3.110 3.215 3.031 3.292 3.426 3.405 3.356 2.993
```


Functions and arguments (suite)
===============================

You can easily create your own function in R. Recommended when you plan to use the same code over and over again.


```r
# Newton-Raphson to find the square root of a number
MySqrt <- function(y) {
  x <- y/2
  while (abs(x*x-y) > 1e-10) {
    x <- (x+y/x)/2
    }
  x
  }
MySqrt(81)
```

```
[1] 9
```

```r
MySqrt(101)
```

```
[1] 10.05
```



Vectorized operation
====================
For loops are notoriously slow in R, and whenever possible, it is preferable to use vectorized operations. Most functions in R are already vectorized.


```r
# Let's generate some uniform [0,10] random numbers
n <- 10000
x  <-  runif(n, 0, 10)
y <- rep(0, n)

library(microbenchmark)
microbenchmark(for(i in 1:n) y[i] <- sqrt(x[i]), sqrt(x), times=10)
```

```
Unit: microseconds
                              expr      min       lq   median       uq
 for (i in 1:n) y[i] <- sqrt(x[i]) 15533.27 15585.87 17164.36 17335.00
                           sqrt(x)    43.64    44.23    45.97    47.03
      max neval
 18099.32    10
    49.01    10
```


The for loop is increadibly slower! 

While you will often hear that R is slow, there are many ways to speed up calculations in R, often by using third party libraries (e.g. data.table, Rcpp). 

The *apply family
===========================================
The `*apply` family of functions is intended to try to solve some of the side effects of `for` loops, such as facilitating it's application to R objects (e.g. lists) and improving efficiency. 


The most common `*apply` functions are 
- `apply`:  Returns a vector or array or list of values obtained by applying a function to margins of an array or matrix.

- `lapply`: apply a function or each element of a list or vector

- `sapply`: a user-friendly version and wrapper of lapply by default returning a vector, matrix or, an array. sapply will try to guess what the output should be based on the input.

For more details have a look at [this](http://www.dummies.com/how-to/content/how-to-use-the-apply-family-of-functions-in-r.html).

Vectorized operation with the *apply family (suite)
===================================================


```r
# Let's generate some uniform [0,10] random numbers
n <- 10000
x  <-  runif(n, 0, 10)
y  <- rep(0, n)
library(microbenchmark)
microbenchmark(for(i in 1:n) y[i] <- MySqrt(x[i]), sapply(x, MySqrt), times=10)
```

```
Unit: milliseconds
                                expr   min    lq median    uq   max neval
 for (i in 1:n) y[i] <- MySqrt(x[i]) 103.5 112.0  119.8 123.8 135.9    10
                   sapply(x, MySqrt) 101.8 104.4  110.0 114.9 122.0    10
```


`*apply` functions are not necessarily faster than `for` loops, but they can be very convenient and usually lead to more compact and more elegant code. 

More efficiency gain can be obtained using compiled code (e.g. C++). R provides multiple ways to call compiled code. In particular, the [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) package can greatly facilitate the use of C++ compiled code.

**Exercise:** Write some code to use the `apply` function on a given matrix.

Reproducibility and literrate programming
=====================
- Approach to programming introduced by Donald Knuth
    - An explanation of a program logic in a plain English, interspersed with chunk of computer code.
- Sweave
    - Create dynamic reports by embedding R code in latex documents
    
knitr vs. Sweave
==================

- Sweave is good but ...
    - Writing latex is painful
    - Output is limited to pdf
- knitr
    - Transparent engine for dynamic report generation
    - knitr allows any input languages (e.g. R, Python and Awk) and any output markup languages
    
knitr: a feature rich package
=============================
- Full control of input, code, and output
    - Fine control over how the code is executed and the ouput is displayed
    - knitr can process input files in various formats: latex, html, R markdown
- [R markdown](http://www.rstudio.com/ide/docs/authoring/using_markdown.html)
    - markdown:  easy-to-read, easy-to-write plain text format that can be converted to html
    - R markdown: markdown + R code chunks
- [R presentation](http://www.rstudio.com/ide/docs/presentations/overview)
- knitr is readily accessible in RStudio

Markdown
========
Markdown is a simple markup language similar to wiki markups

![Markdown in RStudio](http://www.rstudio.com/images/docs/markdownOverview.png)

R Markdown
==========
Mardown with R code chuncks.

![R Mardown](http://www.rstudio.com/images/docs/markdownUntitled.png)

As we've seen, R expressions can also be evaluated inline:

pi=3.1416

knitr and caching
=================

Large data and complex analysis can require significant computing time
    
- Not unusual for an analysis to take a few minutes to an hour, or even more!
- This can result in some performance issues when viewing a report → User frustration
- Why rerun a script when nothing has changed?
- The solution is caching

knitr and caching (suite)
=================

knitr provides powerful caching mechanism:

- cache can be turned on/off for each code chunk
- If caching is on, knitr will check if the code has changed when rerunning a report
- Chunks can be made dependent
- The caching mechanism is flexible can be attached to an R version, an input dataset, a date, etc.

Want to know more about knitr?
=============================
Visit Yihui's webpage: http://yihui.name/knitr/

or buy his book

![Dynamic Documents with R and knitr](http://ecx.images-amazon.com/images/I/41kI1dxXGfL.jpg)

Time for your to work!
=====================
What you need to do:
- Download [R](r-project.org) and [RStudio](rstudio.org)
- Signup for an account on [GitHub](github.com). 
    - Set up your first repository!
- Try [knitr](yihui.name/knitr/) and git within [RStudio](rstudio.org)

We will use RStudio, GitHub and knitr a whole lot throughout this course! 


I expect your to use GitHub/knitr/Rstudio for your homeworks and final project!

