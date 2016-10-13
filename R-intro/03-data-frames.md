---
layout: topic
title: The `data.frame` class
author: Data Carpentry contributors
minutes: 30
---



------------

> ## Learning Objectives
>
> * understand the concept of a `data.frame`
> * use sequences
> * know how to access any element of a `data.frame`

------------

# What are data frames?

`data.frame` is the _de facto_ data structure for most tabular data and what we
use for statistics and plotting.

A `data.frame` is a collection of vectors of identical lengths. Each vector
represents a column, and each vector can be of a different data type (e.g.,
characters, integers, factors). The `str()` function is useful to inspect the
data types of the columns.

A `data.frame` can be created by the functions `read.csv()` or `read.table()`, in
other words, when importing spreadsheets from your hard drive (or the web).

By default, `data.frame` converts (= coerces) columns that contain characters
(i.e., text) into the `factor` data type. Depending on what you want to do with
the data, you may want to keep these columns as `character`. To do so,
`read.csv()` and `read.table()` have an argument called `stringsAsFactors` which
can be set to `FALSE`:


```r
some_data <- read.csv("data/some_file.csv", stringsAsFactors=FALSE)
```

<!--- talk about colClasses argument?, row names?  --->

You can also create `data.frame` manually with the function `data.frame()`. This
function can also take the argument `stringsAsFactors`. Compare the output of
these examples, and compare the difference between when the data are being read
as `character` and when they are being as `factor`.


```r
example_data <- data.frame(animal=c("dog", "cat", "sea cucumber", "sea urchin"),
                           feel=c("furry", "furry", "squishy", "spiny"),
                           weight=c(45, 8, 1.1, 0.8))
str(example_data)
```

```
## 'data.frame':	4 obs. of  3 variables:
##  $ animal: Factor w/ 4 levels "cat","dog","sea cucumber",..: 2 1 3 4
##  $ feel  : Factor w/ 3 levels "furry","spiny",..: 1 1 3 2
##  $ weight: num  45 8 1.1 0.8
```

```r
example_data <- data.frame(animal=c("dog", "cat", "sea cucumber", "sea urchin"),
                           feel=c("furry", "furry", "squishy", "spiny"),
                           weight=c(45, 8, 1.1, 0.8), stringsAsFactors=FALSE)
str(example_data)
```

```
## 'data.frame':	4 obs. of  3 variables:
##  $ animal: chr  "dog" "cat" "sea cucumber" "sea urchin"
##  $ feel  : chr  "furry" "furry" "squishy" "spiny"
##  $ weight: num  45 8 1.1 0.8
```

### Challenge

1. There are a few mistakes in this hand crafted `data.frame`, can you spot and
fix them? Don't hesitate to experiment!


```r
author_book <- data.frame(author_first=c("Charles", "Ernst", "Theodosius"),
                          author_last=c(Darwin, Mayr, Dobzhansky),
                          year=c(1942, 1970))
```

2. Can you predict the class for each of the columns in the following example?


```r
country_climate <- data.frame(country=c("Canada", "Panama", "South Africa", "Australia"),
                              climate=c("cold", "hot", "temperate", "hot/temperate"),
                              temperature=c(10, 30, 18, "15"),
                              northern_hemisphere=c(TRUE, TRUE, FALSE, "FALSE"),
                              has_kangaroo=c(FALSE, FALSE, FALSE, 1))
```

Check your guesses using `str(country_climate)`. Are they what you
expected? Why? Why not?

R coerces (when possible) to the data type that is the least common
denominator and the easiest to coerce to.

# Inspecting `data.frame` objects

We already saw how the functions `head()` and `str()` can be useful to check the
content and the structure of a `data.frame`. Here is a non-exhaustive list of
functions to get a sense of the content/structure of the data.

#### Size:
* `dim()` - returns a vector with the number of rows in the first
  element, and the number of columns as the second element (the
  __dim__ensions of the object)
* `nrow()` - returns the number of rows
* `ncol()` - returns the number of columns

#### Content:
* `head()` - shows the first 6 rows
* `tail()` - shows the last 6 rows

#### Names:
* `names()` - returns the column names (synonym of `colnames()` for `data.frame`
	objects)
* `rownames()` - returns the row names

#### Summary:
* `str()` - structure of the object and information about the class, length and
  content of  each column
* `summary()` - summary statistics for each column

Note: most of these functions are "generic", they can be used on other types of
objects besides `data.frame`.

# Indexing and sequences

If we want to extract one or several values from a vector, we must provide one
or several indices in square brackets, just as we do in math. For instance:


```r
animals <- c("mouse", "rat", "dog", "cat")
animals[2]
```

```
## [1] "rat"
```

```r
animals[c(3, 2)]
```

```
## [1] "dog" "rat"
```

```r
animals[2:4]
```

```
## [1] "rat" "dog" "cat"
```

```r
more_animals <- animals[c(1:3, 2:4)]
more_animals
```

```
## [1] "mouse" "rat"   "dog"   "rat"   "dog"   "cat"
```

R indexes start at 1. Programming languages like Fortran, MATLAB, and R start
counting at 1, because that's what human beings typically do. Languages in the C
family (including C++, Java, Perl, and Python) count from 0 because that's
simpler for computers to do.

`:` is a special function that creates numeric vectors of integer in increasing
or decreasing order, test `1:10` and `10:1` for instance. The function `seq()`
(for __seq__uence) can be used to create more complex patterns:


```r
seq(1, 10, by=2)
```

```
## [1] 1 3 5 7 9
```

```r
seq(5, 10, length.out=3)
```

```
## [1]  5.0  7.5 10.0
```

```r
seq(50, by=5, length.out=10)
```

```
##  [1] 50 55 60 65 70 75 80 85 90 95
```

```r
seq(1, 8, by=3) # sequence stops to stay below upper limit
```

```
## [1] 1 4 7
```

Our survey data frame has rows and columns (it has 2 dimensions), if we want to
extract some specific data from it, we need to specify the "coordinates" we want
from it. Row numbers come first, followed by column numbers.


```r
surveys[1, 1]   # first element in the first column of the data frame
```

```
## [1] 1
```

```r
surveys[1, 6]   # first element in the 6th column
```

```
## [1] NL
## 48 Levels: AB AH AS BA CB CM CQ CS CT CU CV DM DO DS DX NL OL OT OX ... ZL
```

```r
surveys[1:3, 7] # first three elements in the 7th column
```

```
## [1] M M F
## Levels:  F M
```

```r
surveys[3, ]    # the 3rd element for all columns
```

```
##   record_id month day year plot_id species_id sex hindfoot_length weight
## 3         3     7  16 1977       2         DM   F              37     NA
##       genus  species   taxa plot_type
## 3 Dipodomys merriami Rodent   Control
```

```r
surveys[, 8]    # the entire 8th column
```

```
##   [1] 32 33 37 36 35 14 NA 37 34 20 53 38 35 NA 36 36 48 22 NA 48 34 31 36
##  [24] 21 35 31 36 38 NA 52 37 35 36 NA 38 22 35 33 36 36 34 46 36 35 36 35
##  [47] 32 36 17 32 36 26 36 37 36 34 NA 45 33 20 35 35 35 37 34 35 35 32 15
##  [70] 21 36 31 44 12 32 47 NA 16 34 48 14 35 37 35 35 33 11 35 20 35 50 35
##  [93] NA 36 38 36 36 38 37 54
##  [ reached getOption("max.print") -- omitted 34686 entries ]
```

```r
head_surveys <- surveys[1:6, ] # surveys[1:6, ] is equivalent to head(surveys)
```

### Challenge

* The function `nrow()` on a `data.frame` returns the number of
  rows. Use it, in conjuction with `seq()` to create a new
  `data.frame` called `surveys_by_10` that includes every 10th row of
  the survey data frame starting at row 10 (10, 20, 30, ...)



# Indexing by names

For larger datasets, it can be tricky to remember the column number that
corresponds to a particular variable. (Are species names in column 5 or 7? oh,
right... they are in column 6). In some cases, in which column the variable will
be can change if the script you are using adds or removes columns. It's
therefore often better to use column names to refer to a particular variable,
and it makes your code easier to read and your intentions clearer.

You can do operations on a particular column, by selecting it using the `$`
sign. In this case, the entire column is a vector. For instance, to extract all
the weights from our datasets, we can use: `surveys$wgt`. You can use
`names(surveys)` or `colnames(surveys)` to remind yourself of the column names.

In some cases, you may way to select more than one column. You can do this using
the square brackets: `surveys[, c("wgt", "sex")]`.

## Saving and exporting data

We have already seen `saveRDS` and `readRDS` to save and read
serialised R data. These function deal with binary data that can be
used across platforms and R versions.

Note that there is are also two similar functions `save` and
`load`. With the RDS versions, the data is loaded and is then be
assigned into a variable; `load`, on the other hand, loads the
variable in the working space. These functions work with any type of
data (we have seen vectors and data frames so far, but will see more
later).

To export a `data.frame` to text-based spreadsheets, one can use the
`write.csv` function. As it's name implies, it writes the data to a
comma-separated values file.

### Challenge 

Create a new `minisurveys` `data.frame` that contains the 5 first rows
and 3 first columns of `surveys`. 

* Serialise it to disk using either `save` or `saveRDS`, remove it
  from your working space, then load it again (using `load` or
  `readRDS` respectively).

* Export it to `./data/minisurveys.csv`. Delete it and read it back
  in, verifying that it corresponds to the `data.frame` that you
  exported.

# Data frames and tidy data

Data frames can be compared to spreadsheets with additional constraints: 

- All columns must have the same length
- All columns must have the same *type*
- Could can't use colours of font formatting to annotate them
- You can't merge cells

But that's a **advantage**, not a drawback, and here is why: **we want
to know what to expect from the data**, i.e. we want **structured
data** to be able to effectively analyse it (i.e. programme it)
without having to direclty look at it (which anyway becomes
increasingly difficult). We want **tidy data**!

## Tidy data

- Each row contains one observation (one sample)
- Each column documents a single variable
- Each cell contains a single value: never combine multiple pieces of
  information in one cell: split them into different cells
  

![Don't](./img/multiple-info.png)
![Do](./img/single-info.png)

>  Hadley Wickham, Tidy Data, Vol. 59, Issue 10, Sep 2014, Journal of Statistical Software. http://www.jstatsoft.org/v59/i10.

Tidying up data is arguably not the most exciting part of a data
analysis project (as opposed to visualisation, statistical analysis,
...), but it all too often takes a substantial amount of time!

The tidy data rules above are relevant when creating your own
spreadsheets (in Excel, for instance). Assuring the data is tidy and
you avoid information that will be *lost* once the spreadsheet is
exported into a more generic file type (for instance `csv`) will make
it very easily usable into R (or any other data science language).

If you are interested in learning more about good data organisation
and management using spreadsheets, see the
[Data Carpentry Spreadsheet](http://www.datacarpentry.org/spreadsheet-ecology-lesson)
course.

More about [structured data and data structure](./05-datastructures.md). 


