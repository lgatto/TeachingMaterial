---
layout: topic
title: Starting with data
author: Data Carpentry contributors
minutes: 20
---

------------

> ## Learning Objectives
>
> * load external data (CSV files) in memory using the survey table
>  (`surveys.csv`) as an example
> * explore the structure and the content of the data in R
> * understand what are factors and how to manipulate them

------------

# Presentation of the Survey Data

We are studying the species and weight of animals caught in plots in our study
area. The dataset is stored as a `csv` file: each row holds information for a
single animal, and the columns represent:

| Column           | Description                        |
|------------------|------------------------------------|
| record\_id       | Unique id for the observation      |
| month            | month of observation               |
| day              | day of observation                 |
| year             | year of observation                |
| plot\_id         | ID of a particular plot            |
| species\_id      | 2-letter code                      |
| sex              | sex of animal ("M", "F")           |
| hindfoot\_length | length of the hindfoot in mm       |
| weight           | weight of the animal in grams      |
| genus            | genus of animal                    |
| species          | species of animal                  |
| taxa             | e.g. Rodent, Reptile, Bird, Rabbit |
| plot\_type       | type of plot                       |

We are going to use the R function `download.file()` to download the CSV file
that contains the survey data from figshare, and we will use `read.csv()` to
load into memory (as a `data.frame`) the content of the CSV file:

- First, make sure you are in the correct working directory by typing `getwd()`.
- Second, create a new directory within this working directory called `data`. You
can do this by clicking on the new folder icon in RStudio under the file tab, or
by typing `dir.create("data")` at the terminal.
- Third, download the data:


```r
download.file("http://datacarpentry.github.io/dc_zurich/data/portal_data_joined.csv",
              "./data/portal_data_joined.csv")
```

You are now ready to load the data:




```r
surveys <- read.csv('data/portal_data_joined.csv')
```

This statement doesn't produce any output because assignment doesn't display
anything. If we want to check that our data has been loaded, we can print the
variable's value: `surveys`

Alternatively, wrapping an assignment in parentheses will perform the assignment
and display it at the same time.


```r
(surveys <- read.csv('data/portal_data_joined.csv'))
```

Wow... that was a lot of output. At least it means the data loaded
properly. Let's check the top (the first 6 lines) of this `data.frame` using the
function `head()`:


```r
head(surveys)
```

```
##   record_id month day year plot_id species_id sex hindfoot_length weight
## 1         1     7  16 1977       2         NL   M              32     NA
## 2         2     7  16 1977       3         NL   M              33     NA
## 3         3     7  16 1977       2         DM   F              37     NA
## 4         4     7  16 1977       7         DM   M              36     NA
## 5         5     7  16 1977       3         DM   M              35     NA
## 6         6     7  16 1977       1         PF   M              14     NA
##         genus  species   taxa                plot_type
## 1     Neotoma albigula Rodent                  Control
## 2     Neotoma albigula Rodent Long-term Krat Exclosure
## 3   Dipodomys merriami Rodent                  Control
## 4   Dipodomys merriami Rodent         Rodent Exclosure
## 5   Dipodomys merriami Rodent Long-term Krat Exclosure
## 6 Perognathus   flavus Rodent        Spectab exclosure
```

Let's now check the __str__ucture of this `data.frame` in more details with the
function `str()`:


```r
str(surveys)
```

```
## 'data.frame':	34786 obs. of  13 variables:
##  $ record_id      : int  1 2 3 4 5 6 7 8 9 10 ...
##  $ month          : int  7 7 7 7 7 7 7 7 7 7 ...
##  $ day            : int  16 16 16 16 16 16 16 16 16 16 ...
##  $ year           : int  1977 1977 1977 1977 1977 1977 1977 1977 1977 1977 ...
##  $ plot_id        : int  2 3 2 7 3 1 2 1 1 6 ...
##  $ species_id     : Factor w/ 48 levels "AB","AH","AS",..: 16 16 12 12 12 23 22 12 12 23 ...
##  $ sex            : Factor w/ 3 levels "","F","M": 3 3 2 3 3 3 2 3 2 2 ...
##  $ hindfoot_length: int  32 33 37 36 35 14 NA 37 34 20 ...
##  $ weight         : int  NA NA NA NA NA NA NA NA NA NA ...
##  $ genus          : Factor w/ 26 levels "Ammodramus","Ammospermophilus",..: 13 13 11 11 11 15 16 11 11 15 ...
##  $ species        : Factor w/ 40 levels "albigula","audubonii",..: 1 1 23 23 23 9 8 23 23 9 ...
##  $ taxa           : Factor w/ 4 levels "Bird","Rabbit",..: 4 4 4 4 4 4 4 4 4 4 ...
##  $ plot_type      : Factor w/ 5 levels "Control","Long-term Krat Exclosure",..: 1 2 1 3 2 5 1 5 5 4 ...
```

### Challenge

Based on the output of `str(surveys)`, can you answer the following questions?

* What is the class of the object `surveys`?
* How many rows and how many columns are in this object?
* How many species have been recorded during these surveys?



As you can see, the columns `species` and `sex` are of a special class called
`factor`. Before we learn more about the `data.frame` class, we are going to
talk about factors. They are very useful but not necessarily intuitive, and
therefore require some attention.


## Factors

Factors are used to represent categorical data. Factors can be ordered or
unordered and are an important class for statistical analysis and for plotting.

Factors are stored as integers, and have labels associated with these unique
integers. While factors look (and often behave) like character vectors, they are
actually integers under the hood, and you need to be careful when treating them
like strings.

Once created, factors can only contain a pre-defined set values, known as
*levels*. By default, R always sorts *levels* in alphabetical order. For
instance, if you have a factor with 2 levels:


```r
sex <- factor(c("male", "female", "female", "male"))
```

R will assign `1` to the level `"female"` and `2` to the level `"male"` (because
`f` comes before `m`, even though the first element in this vector is
`"male"`). You can check this by using the function `levels()`, and check the
number of levels using `nlevels()`:


```r
levels(sex)
```

```
## [1] "female" "male"
```

```r
nlevels(sex)
```

```
## [1] 2
```

Sometimes, the order of the factors does not matter, other times you might want
to specify the order because it is meaningful (e.g., "low", "medium", "high") or
it is required by particular type of analysis. Additionally, specifying the
order of the levels allows to compare levels:


```r
food <- factor(c("low", "high", "medium", "high", "low", "medium", "high"))
levels(food)
```

```
## [1] "high"   "low"    "medium"
```

```r
food <- factor(food, levels=c("low", "medium", "high"))
levels(food)
```

```
## [1] "low"    "medium" "high"
```

```r
min(food) ## doesn't work
```

```
## Error in Summary.factor(structure(c(1L, 3L, 2L, 3L, 1L, 2L, 3L), .Label = c("low", : 'min' not meaningful for factors
```

```r
food <- factor(food, levels=c("low", "medium", "high"), ordered=TRUE)
levels(food)
```

```
## [1] "low"    "medium" "high"
```

```r
min(food) ## works!
```

```
## [1] low
## Levels: low < medium < high
```

In R's memory, these factors are represented by numbers (1, 2, 3). They are
better than using simple integer labels because factors are self describing:
`"low"`, `"medium"`, and `"high"`" is more descriptive than `1`, `2`, `3`. Which
is low?  You wouldn't be able to tell with just integer data. Factors have this
information built in. It is particularly helpful when there are many levels
(like the species in our example data set).

### Converting factors

If you need to convert a factor to a character vector, simply use
`as.character(x)`.

Converting a factor to a numeric vector is however a little trickier, and you
have to go via a character vector. Compare:


```r
f <- factor(c(1, 5, 10, 2))
as.numeric(f)               ## wrong! and there is no warning...
```

```
## [1] 1 3 4 2
```

```r
as.numeric(as.character(f)) ## works...
```

```
## [1]  1  5 10  2
```

```r
as.numeric(levels(f))[f]    ## The recommended way.
```

```
## [1]  1  5 10  2
```

### Challenge

The function `table()` tabulates observations and can be used to
create bar plots quickly. For instance, how can you recreate this plot
but by having `control` being listed last instead of first?


```r
exprmt <- factor(c("treat1", "treat2", "treat1", "treat3", "treat1", "control",
                   "treat1", "treat2", "treat3"))
table(exprmt)
```

```
## exprmt
## control  treat1  treat2  treat3 
##       1       4       2       2
```

```r
barplot(table(exprmt))
```

![plot of chunk wrong-order](figure/wrong-order-1.png)


