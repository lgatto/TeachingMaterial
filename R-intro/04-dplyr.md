---
layout: topic
title: Aggregating and analyzing data with dplyr
author: Data Carpentry contributors
---



------------

# Data manipulation using `dplyr`

Bracket subsetting is handy, but it can be cumbersome and difficult to read,
especially for complicated operations. Enter `dplyr`. `dplyr` is a package for
making data manipulation easier.

Packages in R are basically sets of additional functions that let you do more
stuff in R. The functions we've been using, like `str()`, come built into R;
packages give you access to more functions. You need to install a package and
then load it to be able to use it.


```r
install.packages("dplyr") ## install
```

You might get asked to choose a **CRAN mirror** -- this is basically
asking you to choose a site to download the package from. The choice
doesn't matter too much. A local CRAN mirror is a good choice.



```r
library("dplyr")          ## load
```

You only need to install a package once per computer, but you need to load it
every time you open a new R session and want to use that package.

## What is `dplyr`?

> The package `dplyr` provides **easy tools for the most common data
> manipulation tasks**. It is built to work directly with **tidy**
> data frames: `dplyr`'s function take tidy data, transform it, and
> return a new tidy data.

The thinking behind it was largely inspired by the package `plyr`
which has been in use for some time but suffered from being slow in
some cases.` dplyr` addresses this by porting much of the computation
to C++. An additional feature is the ability to work with data stored
directly in an external database. The benefits of doing this are that
the data can be managed natively in a relational database, queries can
be conducted on that database, and only the results of the query
returned.

This addresses a common problem with R in that all operations are conducted in
memory and thus the amount of data you can work with is limited by available
memory. The database connections essentially remove that limitation in that you
can have a database of many 100s GB, conduct queries on it directly and pull
back just what you need for analysis in R.

> We're going to learn some of the most common `dplyr` functions:
> `select()`, `filter()`, `mutate()`, `group_by()`, and `summarize()`
> (or `summarise()`).

### Selecting columns and filtering rows

To select columns of a data frame, use `select()`. The first argument
to this function is the data frame (`surveys`), and the subsequent
arguments are the columns to keep.


```r
select(surveys, plot_id, species_id, weight)
```

To choose rows, use `filter()`:


```r
filter(surveys, year == 1995)
```

```
##      record_id month day year plot_id species_id sex hindfoot_length
## 1        21993     1  11 1995      18         PF   F              16
## 2        21994     1  11 1995      12         DO   M              36
## 3        21995     1  11 1995       2         DO   M              36
## 4        21996     1  11 1995      21         PF   F              14
## 5        21997     1  11 1995      24         RM   M              15
## 6        21998     1  11 1995       1         DM   M              38
## 7        21999     1  11 1995      19         PF   F              15
##      weight            genus         species    taxa
## 1         7      Perognathus          flavus  Rodent
## 2        47        Dipodomys           ordii  Rodent
## 3        51        Dipodomys           ordii  Rodent
## 4         7      Perognathus          flavus  Rodent
## 5        10  Reithrodontomys       megalotis  Rodent
## 6        46        Dipodomys        merriami  Rodent
## 7         8      Perognathus          flavus  Rodent
##                      plot_type
## 1    Short-term Krat Exclosure
## 2                      Control
## 3                      Control
## 4     Long-term Krat Exclosure
## 5             Rodent Exclosure
## 6            Spectab exclosure
## 7     Long-term Krat Exclosure
##  [ reached getOption("max.print") -- omitted 1173 rows ]
```

### Pipes

**But what if you wanted to select and filter?** There are three ways to do this:

- Use intermediate steps, nested functions, or pipes. With the
  intermediate steps, you essentially create a temporary data frame
  and use that as input to the next function. This can clutter up your
  workspace with lots of objects.

- You can also nest functions (i.e. one function inside of
  another). This is handy, but can be difficult to read if too many
  functions are nested as the process from inside out.

- The last option, pipes, are a fairly recent addition to R.

Pipes let you take the output of one function and send it directly to
the next, which is useful when you need to many things to the same
data set.  Pipes in R look like `%>%` and are made available via the
`magrittr` package installed as part of `dplyr`.


```r
surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight)
```

```
##    species_id sex weight
## 1          PF   M      4
## 2          PF   F      4
## 3          PF          4
## 4          PF   F      4
## 5          PF   F      4
## 6          RM   M      4
## 7          RM   F      4
## 8          RM   M      4
## 9          RM   M      4
## 10         RM   M      4
## 11         RM   M      4
## 12         RM   F      4
## 13         RM   M      4
## 14         RM   M      4
## 15         RM   M      4
## 16         PF   M      4
## 17         PP   M      4
```

In the above we use the pipe to send the `surveys` data set first
through `filter`, to keep rows where `weight` was less than 5, and
then through `select` to keep the `species`, `sex` and `weight`'
columns. When the data frame is being passed to the `filter()` and
`select()` functions through a pipe, we don't need to include it as an
argument to these functions anymore.

If we wanted to create a new object with this smaller version of the
data we could do so by assigning it a new name:


```r
surveys_sml <- surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight)

surveys_sml
```

```
##    species_id sex weight
## 1          PF   M      4
## 2          PF   F      4
## 3          PF          4
## 4          PF   F      4
## 5          PF   F      4
## 6          RM   M      4
## 7          RM   F      4
## 8          RM   M      4
## 9          RM   M      4
## 10         RM   M      4
## 11         RM   M      4
## 12         RM   F      4
## 13         RM   M      4
## 14         RM   M      4
## 15         RM   M      4
## 16         PF   M      4
## 17         PP   M      4
```

### Challenge

> Using pipes, subset the data to include rows before 1995. Retain
> columns `year`, `sex`, and `weight.`





### Mutate

Frequently you'll want to create new columns based on the values in existing
columns, for example to do unit conversions or find the ratio of values in two
columns. For this we'll use `mutate()`.

To create a new column of weight in kg:


```r
surveys %>%
  mutate(weight_kg = weight / 1000)
```

```
##       record_id month day year plot_id species_id sex hindfoot_length
## 1             1     7  16 1977       2         NL   M              32
## 2             2     7  16 1977       3         NL   M              33
## 3             3     7  16 1977       2         DM   F              37
## 4             4     7  16 1977       7         DM   M              36
## 5             5     7  16 1977       3         DM   M              35
## 6             6     7  16 1977       1         PF   M              14
## 7             7     7  16 1977       2         PE   F              NA
##       weight            genus         species    taxa
## 1         NA          Neotoma        albigula  Rodent
## 2         NA          Neotoma        albigula  Rodent
## 3         NA        Dipodomys        merriami  Rodent
## 4         NA        Dipodomys        merriami  Rodent
## 5         NA        Dipodomys        merriami  Rodent
## 6         NA      Perognathus          flavus  Rodent
## 7         NA       Peromyscus        eremicus  Rodent
##                       plot_type weight_kg
## 1                       Control        NA
## 2      Long-term Krat Exclosure        NA
## 3                       Control        NA
## 4              Rodent Exclosure        NA
## 5      Long-term Krat Exclosure        NA
## 6             Spectab exclosure        NA
## 7                       Control        NA
##  [ reached getOption("max.print") -- omitted 34779 rows ]
```

If this runs off your screen and you just want to see the first few rows, you
can use a pipe to view the `head()` of the data (pipes work with non-dplyr
functions too, as long as the `dplyr` or `magrittr` packages are loaded).


```r
surveys %>%
  mutate(weight_kg = weight / 1000) %>%
  head
```

```
##   record_id month day year plot_id species_id sex hindfoot_length weight
## 1         1     7  16 1977       2         NL   M              32     NA
## 2         2     7  16 1977       3         NL   M              33     NA
## 3         3     7  16 1977       2         DM   F              37     NA
## 4         4     7  16 1977       7         DM   M              36     NA
## 5         5     7  16 1977       3         DM   M              35     NA
## 6         6     7  16 1977       1         PF   M              14     NA
##         genus  species   taxa                plot_type weight_kg
## 1     Neotoma albigula Rodent                  Control        NA
## 2     Neotoma albigula Rodent Long-term Krat Exclosure        NA
## 3   Dipodomys merriami Rodent                  Control        NA
## 4   Dipodomys merriami Rodent         Rodent Exclosure        NA
## 5   Dipodomys merriami Rodent Long-term Krat Exclosure        NA
## 6 Perognathus   flavus Rodent        Spectab exclosure        NA
```

The first few rows are full of NAs, so if we wanted to remove those we could
insert a `filter()` in this chain:


```r
surveys %>%
  mutate(weight_kg = weight / 1000) %>%
  filter(!is.na(weight)) %>%
  head
```

```
##   record_id month day year plot_id species_id sex hindfoot_length weight
## 1        63     8  19 1977       3         DM   M              35     40
## 2        64     8  19 1977       7         DM   M              37     48
## 3        65     8  19 1977       4         DM   F              34     29
## 4        66     8  19 1977       4         DM   F              35     46
## 5        67     8  19 1977       7         DM   M              35     36
## 6        68     8  19 1977       8         DO   F              32     52
##       genus  species   taxa                plot_type weight_kg
## 1 Dipodomys merriami Rodent Long-term Krat Exclosure     0.040
## 2 Dipodomys merriami Rodent         Rodent Exclosure     0.048
## 3 Dipodomys merriami Rodent                  Control     0.029
## 4 Dipodomys merriami Rodent                  Control     0.046
## 5 Dipodomys merriami Rodent         Rodent Exclosure     0.036
## 6 Dipodomys    ordii Rodent                  Control     0.052
```

`is.na()` is a function that determines whether something is or is not an `NA`.
The `!` symbol negates it, so we're asking for everything that is not an `NA`.

### Split-apply-combine data analysis and the summarize() function

Many data analysis tasks can be approached using the "split-apply-combine"
paradigm:

1. split the data into groups
2. apply some analysis to each group, and
3. then combine the results.

`dplyr` makes this very easy through the use of the `group_by()`
function. `group_by()` splits the data into groups upon which some
operations can be run. For example, if we wanted to group by sex and
find the number of rows of data for each sex, we would do:


```r
surveys %>%
  group_by(sex) %>%
  tally()
```

```
## # A tibble: 3 × 2
##      sex     n
##   <fctr> <int>
## 1         1748
## 2      F 15690
## 3      M 17348
```

`group_by()` is often used together with `summarize()` which collapses each
group into a single-row summary of that group. So to view mean `weight` by sex:


```r
surveys %>%
  group_by(sex) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))
```

```
## # A tibble: 3 × 2
##      sex mean_weight
##   <fctr>       <dbl>
## 1           64.74257
## 2      F    42.17055
## 3      M    42.99538
```

You can group by multiple columns too:


```r
surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))
```

```
## Source: local data frame [92 x 3]
## Groups: sex [?]
## 
##       sex species_id mean_weight
##    <fctr>     <fctr>       <dbl>
## 1                 AB         NaN
## 2                 AH         NaN
## 3                 AS         NaN
## 4                 BA         NaN
## 5                 CB         NaN
## 6                 CM         NaN
## 7                 CQ         NaN
## 8                 CS         NaN
## 9                 CT         NaN
## 10                CU         NaN
## # ... with 82 more rows
```

Looks like most of these species were never weighed. We could then discard rows
where `mean_weight` is `NaN` (`NaN` refers to "Not a Number") using `filter()`:


```r
surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE)) %>%
  filter(!is.nan(mean_weight))
```

```
## Source: local data frame [64 x 3]
## Groups: sex [3]
## 
##       sex species_id mean_weight
##    <fctr>     <fctr>       <dbl>
## 1                 DM    38.28571
## 2                 DO    50.66667
## 3                 DS   120.00000
## 4                 NL   167.68750
## 5                 OL    29.00000
## 6                 OT    21.20000
## 7                 PB    30.60000
## 8                 PE    17.66667
## 9                 PF     6.00000
## 10                PI    18.00000
## # ... with 54 more rows
```

All of a sudden this isn't running of the screen anymore. That's because `dplyr`
has changed our `data.frame` to a `tbl_df`. This is a data structure that's very
similar to a data frame; for our purposes the only difference is that it won't
automatically show tons of data going off the screen.

You can also summarize multiple variables at the same time:


```r
surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE),
            min_weight = min(weight, na.rm = TRUE)) %>%
  filter(!is.nan(mean_weight))
```

```
## Source: local data frame [64 x 4]
## Groups: sex [3]
## 
##       sex species_id mean_weight min_weight
##    <fctr>     <fctr>       <dbl>      <int>
## 1                 DM    38.28571         24
## 2                 DO    50.66667         44
## 3                 DS   120.00000         78
## 4                 NL   167.68750         83
## 5                 OL    29.00000         21
## 6                 OT    21.20000         18
## 7                 PB    30.60000         20
## 8                 PE    17.66667         17
## 9                 PF     6.00000          4
## 10                PI    18.00000         18
## # ... with 54 more rows
```

### Challenge

> How many times was each `plot_type` surveyed?




> Use `group_by()` and `summarize()` to find the mean, min, and max hindfoot
> length for each species.



> What was the heaviest animal measured in each year? Return the
> columns `year`, `genus`, `species`, and `weight`. See also the
> cheatsheet below for help.



### Summary

* `select`: select columns
* `filter`: select rows
* `mutate`: create new columns
* `group_by`: split-apply-combine
* `summarise`: collapse each group into a single-row summary of that
  group
* `%>%`: to pipe (chain) operations

[Handy dplyr cheatsheet](http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)

*Much of this lesson was copied or adapted from Jeff Hollister's [materials](http://usepa.github.io/introR/2015/01/14/03-Clean/)*
