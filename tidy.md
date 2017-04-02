# Tidyverse

The **tidyverse**, popularised by Hadley Wickham, refers to the
utilisation of **tidy data** and **tools that preserve tidyness**.

We are going to 

* Learn how tidy data is defined for dataframes
* Learn some tidyverse tools
* Extend the notion of tidy data to Bioconductor's rich semantic

From http://tidyverse.org/:

> The tidyverse is a set of packages that work in harmony because they
> share common data representations and API design. The tidyverse
> package is designed to make it easy to install and load core
> packages from the tidyverse in a single command.



```r
library("tidyverse")
```

```
## Loading tidyverse: ggplot2
## Loading tidyverse: tibble
## Loading tidyverse: tidyr
## Loading tidyverse: readr
## Loading tidyverse: purrr
## Loading tidyverse: dplyr
```

```
## Conflicts with tidy packages ----------------------------------------------
```

```
## filter(): dplyr, stats
## lag():    dplyr, stats
```

(but we are only going to use a subet of these packages)

## Tidy data

Just as

> Happy families are all alike; every unhappy family is unhappy in its
> own way. – Leo Tolstoy

when it comes to data

> Tidy datasets are all alike, but every messy dataset is messy in its
> own way. – Hadley Wickham

The reason why tidy data is important is that we waste a lot of time
in cleaning *messy* data, i.e. tiding it up to get it in a format that
leads to data analysis and visualisation.

One reason why data is often messy is that its structure is meant to
make collection of the data easy.

The standard representation of data is in the form of a table. Tidy
tables are tables where:

1. Each variable is in a column.
2. Each observation is a row.
3. Each value is a cell.

### Examples

Beyond messy badly formatted data (from the
[Data Carpentry *Spreadsheet lessons*](http://www.datacarpentry.org/spreadsheet-ecology-lesson/))

![badly formatted data](./figs/multiple-info.png)


Untidy and tidyfied data (from [Wickham](http://www.jstatsoft.org/v59/i10) 2014.)

![untidy data](./figs/tidy-ex-1a.png)![tidy data](./figs/tidy-ex-1b.png)

From the `tidyr` package: country, year, cases and population from the
World Health Organization Global Tuberculosis Report, organised in
four different ways.


```r
table1
```

```
## # A tibble: 6 × 4
##       country  year  cases population
##         <chr> <int>  <int>      <int>
## 1 Afghanistan  1999    745   19987071
## 2 Afghanistan  2000   2666   20595360
## 3      Brazil  1999  37737  172006362
## 4      Brazil  2000  80488  174504898
## 5       China  1999 212258 1272915272
## 6       China  2000 213766 1280428583
```

```r
table2
```

```
## # A tibble: 12 × 4
##        country  year       type      count
##          <chr> <int>      <chr>      <int>
## 1  Afghanistan  1999      cases        745
## 2  Afghanistan  1999 population   19987071
## 3  Afghanistan  2000      cases       2666
## 4  Afghanistan  2000 population   20595360
## 5       Brazil  1999      cases      37737
## 6       Brazil  1999 population  172006362
## 7       Brazil  2000      cases      80488
## 8       Brazil  2000 population  174504898
## 9        China  1999      cases     212258
## 10       China  1999 population 1272915272
## 11       China  2000      cases     213766
## 12       China  2000 population 1280428583
```

```r
table3
```

```
## # A tibble: 6 × 3
##       country  year              rate
## *       <chr> <int>             <chr>
## 1 Afghanistan  1999      745/19987071
## 2 Afghanistan  2000     2666/20595360
## 3      Brazil  1999   37737/172006362
## 4      Brazil  2000   80488/174504898
## 5       China  1999 212258/1272915272
## 6       China  2000 213766/1280428583
```

```r
table4a
```

```
## # A tibble: 3 × 3
##       country `1999` `2000`
## *       <chr>  <int>  <int>
## 1 Afghanistan    745   2666
## 2      Brazil  37737  80488
## 3       China 212258 213766
```

```r
table4b
```

```
## # A tibble: 3 × 3
##       country     `1999`     `2000`
## *       <chr>      <int>      <int>
## 1 Afghanistan   19987071   20595360
## 2      Brazil  172006362  174504898
## 3       China 1272915272 1280428583
```

## Tidy tools: manipulating and analyzing data with `dplyr`

Credit: This material is based on the Data Carpentry
[*R for data analysis and visualization of Ecological Data* material](http://www.datacarpentry.org/R-ecology-lesson/index.html)

### The data

The *survey data* provides the species and weight of animals caught in
plots in various study area. The dataset is stored as a comma
separated value (CSV) file.  Each row holds information for a single
animal, and the columns represent:

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

It is available online
[https://ndownloader.figshare.com/files/2292169](https://ndownloader.figshare.com/files/2292169). We read it directly from figshare using `dplyr::read_csv`.


```r
surveys <- read_csv("https://ndownloader.figshare.com/files/2292169")
```

```
## Parsed with column specification:
## cols(
##   record_id = col_integer(),
##   month = col_integer(),
##   day = col_integer(),
##   year = col_integer(),
##   plot_id = col_integer(),
##   species_id = col_character(),
##   sex = col_character(),
##   hindfoot_length = col_integer(),
##   weight = col_integer(),
##   genus = col_character(),
##   species = col_character(),
##   taxa = col_character(),
##   plot_type = col_character()
## )
```

```r
head(surveys)
```

```
## # A tibble: 6 × 13
##   record_id month   day  year plot_id species_id   sex hindfoot_length
##       <int> <int> <int> <int>   <int>      <chr> <chr>           <int>
## 1         1     7    16  1977       2         NL     M              32
## 2        72     8    19  1977       2         NL     M              31
## 3       224     9    13  1977       2         NL  <NA>              NA
## 4       266    10    16  1977       2         NL  <NA>              NA
## 5       349    11    12  1977       2         NL  <NA>              NA
## 6       363    11    12  1977       2         NL  <NA>              NA
## # ... with 5 more variables: weight <int>, genus <chr>, species <chr>,
## #   taxa <chr>, plot_type <chr>
```

### Piping

When running complex operations on data, we often have to either
create temporary variables, or nest functions. An alternative is
piping operations using the `%>%` operator from the `magrittr`
package.


```r
library("magrittr")
dim(surveys)
```

```
## [1] 34786    13
```

```r
surveys %>% dim
```

```
## [1] 34786    13
```

This is particularly suited to tidy data and tidy tools. Instead of 

```
tidy_data_new <- tidy_tool1(tidy_data)
tidy_data_new <- tidy_tool2(tidy_data_new)
new_tidy_data_new <- tidy_tool3(tidy_data_new)
```
we have
```
tidy_data_new <- tidy_data %>% tidy_tool1 %>% tidy_tool2 %>% tidy_tools3
```

### Selecting variables (columns) with `select`


```r
surveys %>% select(species, plot_id, species_id)
```

```
## # A tibble: 34,786 × 3
##     species plot_id species_id
##       <chr>   <int>      <chr>
## 1  albigula       2         NL
## 2  albigula       2         NL
## 3  albigula       2         NL
## 4  albigula       2         NL
## 5  albigula       2         NL
## 6  albigula       2         NL
## 7  albigula       2         NL
## 8  albigula       2         NL
## 9  albigula       2         NL
## 10 albigula       2         NL
## # ... with 34,776 more rows
```

### Selecting observations (rows) with `filter`


```r
surveys %>% filter(year == 1995)
```

```
## # A tibble: 1,180 × 13
##    record_id month   day  year plot_id species_id   sex hindfoot_length
##        <int> <int> <int> <int>   <int>      <chr> <chr>           <int>
## 1      22314     6     7  1995       2         NL     M              34
## 2      22728     9    23  1995       2         NL     F              32
## 3      22899    10    28  1995       2         NL     F              32
## 4      23032    12     2  1995       2         NL     F              33
## 5      22003     1    11  1995       2         DM     M              37
## 6      22042     2     4  1995       2         DM     F              36
## 7      22044     2     4  1995       2         DM     M              37
## 8      22105     3     4  1995       2         DM     F              37
## 9      22109     3     4  1995       2         DM     M              37
## 10     22168     4     1  1995       2         DM     M              36
## # ... with 1,170 more rows, and 5 more variables: weight <int>,
## #   genus <chr>, species <chr>, taxa <chr>, plot_type <chr>
```

Pipes come handy when we want to `select` and `filter`:


```r
surveys %>%
  filter(weight < 5) %>%
  select(species, sex, weight)
```

```
## # A tibble: 17 × 3
##         species   sex weight
##           <chr> <chr>  <int>
## 1        flavus     F      4
## 2        flavus     F      4
## 3        flavus     M      4
## 4     megalotis     F      4
## 5     megalotis     M      4
## 6        flavus  <NA>      4
## 7  penicillatus     M      4
## 8     megalotis     M      4
## 9     megalotis     M      4
## 10    megalotis     M      4
## 11       flavus     M      4
## 12       flavus     F      4
## 13    megalotis     M      4
## 14    megalotis     M      4
## 15    megalotis     F      4
## 16    megalotis     M      4
## 17    megalotis     M      4
```

And to save the final, left-most, result


```r
sml <- surveys %>%
    filter(weight < 5) %>%
    select(species, sex, weight)
```

> ### Challenge 
>
>  Using pipes, subset the `survey` data to include individuals collected before
>  1995 and retain only the columns `year`, `sex`, and `weight`.

<!---

```r
## Answer
surveys %>%
    filter(year < 1995) %>%
    select(year, sex, weight)
```
--->

### Adding variables with `mutate`


```r
surveys %>%
    mutate(weight_kg = weight / 1000) %>%
    select(species, year, weight, weight_kg) 
```

```
## # A tibble: 34,786 × 4
##     species  year weight weight_kg
##       <chr> <int>  <int>     <dbl>
## 1  albigula  1977     NA        NA
## 2  albigula  1977     NA        NA
## 3  albigula  1977     NA        NA
## 4  albigula  1977     NA        NA
## 5  albigula  1977     NA        NA
## 6  albigula  1977     NA        NA
## 7  albigula  1977     NA        NA
## 8  albigula  1978     NA        NA
## 9  albigula  1978    218     0.218
## 10 albigula  1978     NA        NA
## # ... with 34,776 more rows
```



```r
surveys %>%
  filter(!is.na(weight)) %>%
    mutate(weight_kg = weight / 1000) %>%
    select(species, year, weight, weight_kg) 
```

```
## # A tibble: 32,283 × 4
##     species  year weight weight_kg
##       <chr> <int>  <int>     <dbl>
## 1  albigula  1978    218     0.218
## 2  albigula  1978    204     0.204
## 3  albigula  1978    200     0.200
## 4  albigula  1978    199     0.199
## 5  albigula  1978    197     0.197
## 6  albigula  1978    218     0.218
## 7  albigula  1979    166     0.166
## 8  albigula  1979    184     0.184
## 9  albigula  1979    206     0.206
## 10 albigula  1979    274     0.274
## # ... with 32,273 more rows
```

> ### Challenge
>
>  Create a new data frame from the `survey` data that meets the following
>  criteria: contains only the `species_id` column and a new column called
>  `hindfoot_half` containing values that are half the `hindfoot_length` values.
>  In this `hindfoot_half` column, there are no `NA`s and all values are less
>  than 30.


<!---

```r
## Answer
surveys_hindfoot_half <- surveys %>%
    filter(!is.na(hindfoot_length)) %>%
    mutate(hindfoot_half = hindfoot_length / 2) %>%
    filter(hindfoot_half < 30) %>%
    select(species_id, hindfoot_half)
```
--->

### Split-apply-combine data analysis and the summarize() function

Many data analysis tasks can be approached using the *split-apply-combine*
paradigm: split the data into groups, apply some analysis to each group, and
then combine the results. `dplyr` makes this very easy through the use of the
`group_by()` function.


#### The `summarize()` function

`group_by()` is often used together with `summarize()`, which collapses each
group into a single-row summary of that group.  `group_by()` takes as arguments
the column names that contain the **categorical** variables for which you want
to calculate the summary statistics. So to view the mean `weight` by sex:


```r
surveys %>%
  group_by(sex) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))
```

```
## # A tibble: 3 × 2
##     sex mean_weight
##   <chr>       <dbl>
## 1     F    42.17055
## 2     M    42.99538
## 3  <NA>    64.74257
```

You can also group by multiple columns:


```r
surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))
```

```
## Source: local data frame [92 x 3]
## Groups: sex [?]
## 
##      sex species_id mean_weight
##    <chr>      <chr>       <dbl>
## 1      F         BA     9.16129
## 2      F         DM    41.60968
## 3      F         DO    48.53125
## 4      F         DS   117.74955
## 5      F         NL   154.28221
## 6      F         OL    31.06582
## 7      F         OT    24.83090
## 8      F         OX    21.00000
## 9      F         PB    30.21088
## 10     F         PE    22.82218
## # ... with 82 more rows
```

When grouping both by `sex` and `species_id`, the first rows are for individuals
that escaped before their sex could be determined and weighted. You may notice
that the last column does not contain `NA` but `NaN` (which refers to "Not a
Number"). To avoid this, we can remove the missing values for weight before we
attempt to calculate the summary statistics on weight. Because the missing
values are removed, we can omit `na.rm = TRUE` when computing the mean:


```r
surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight))
```

```
## Source: local data frame [64 x 3]
## Groups: sex [?]
## 
##      sex species_id mean_weight
##    <chr>      <chr>       <dbl>
## 1      F         BA     9.16129
## 2      F         DM    41.60968
## 3      F         DO    48.53125
## 4      F         DS   117.74955
## 5      F         NL   154.28221
## 6      F         OL    31.06582
## 7      F         OT    24.83090
## 8      F         OX    21.00000
## 9      F         PB    30.21088
## 10     F         PE    22.82218
## # ... with 54 more rows
```

You may also have noticed that the output from these calls doesn't run off the
screen anymore. That's because `dplyr` has changed our `data.frame` to a
`tbl_df`. The `tbl` data structure is very similar to a data frame; for our
purposes the only difference is that, in addition to displaying the data type
of each column under its name, it only prints the first few rows of data and
only as many columns as fit on one screen. If you want to display more data, you
use the `print()` function at the end of your chain with the argument `n`
specifying the number of rows to display:


```r
surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight)) %>%
  print(n = 15)
```

```
## Source: local data frame [64 x 3]
## Groups: sex [?]
## 
##      sex species_id mean_weight
##    <chr>      <chr>       <dbl>
## 1      F         BA    9.161290
## 2      F         DM   41.609685
## 3      F         DO   48.531250
## 4      F         DS  117.749548
## 5      F         NL  154.282209
## 6      F         OL   31.065817
## 7      F         OT   24.830904
## 8      F         OX   21.000000
## 9      F         PB   30.210884
## 10     F         PE   22.822183
## 11     F         PF    7.974394
## 12     F         PH   30.850000
## 13     F         PL   19.312500
## 14     F         PM   22.125668
## 15     F         PP   17.180670
## # ... with 49 more rows
```

Once the data are grouped, you can also summarize multiple variables at the same
time (and not necessarily on the same variable). For instance, we could add a
column indicating the minimum weight for each species for each sex:


```r
surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight),
            min_weight = min(weight))
```

```
## Source: local data frame [64 x 4]
## Groups: sex [?]
## 
##      sex species_id mean_weight min_weight
##    <chr>      <chr>       <dbl>      <int>
## 1      F         BA     9.16129          6
## 2      F         DM    41.60968         10
## 3      F         DO    48.53125         12
## 4      F         DS   117.74955         45
## 5      F         NL   154.28221         32
## 6      F         OL    31.06582         10
## 7      F         OT    24.83090          5
## 8      F         OX    21.00000         20
## 9      F         PB    30.21088         12
## 10     F         PE    22.82218         11
## # ... with 54 more rows
```


#### Tallying

When working with data, it is also common to want to know the number of
observations found for each factor or combination of factors. For this, `dplyr`
provides `tally()`. For example, if we wanted to group by sex and find the
number of rows of data for each sex, we would do:


```r
surveys %>%
  group_by(sex) %>%
  tally
```

```
## # A tibble: 3 × 2
##     sex     n
##   <chr> <int>
## 1     F 15690
## 2     M 17348
## 3  <NA>  1748
```

Here, `tally()` is the action applied to the groups created by `group_by()` and
counts the total number of records for each category.

> ### Challenge
>
> 1. How many individuals were caught in each `plot_type` surveyed?
>
> 2. Use `group_by()` and `summarize()` to find the mean, min, and max hindfoot
> length for each species (using `species_id`).
>
> 3. What was the heaviest animal measured in each year? Return the columns `year`,
> `genus`, `species_id`, and `weight`.
>
> 4. You saw above how to count the number of individuals of each `sex` using a
> combination of `group_by()` and `tally()`. How could you get the same result
> using `group_by()` and `summarize()`? Hint: see `?n`.


<!---

```
## # A tibble: 5 × 2
##                   plot_type     n
##                       <chr> <int>
## 1                   Control 15611
## 2  Long-term Krat Exclosure  5118
## 3          Rodent Exclosure  4233
## 4 Short-term Krat Exclosure  5906
## 5         Spectab exclosure  3918
```

```
## # A tibble: 25 × 4
##    species_id mean_hindfoot_length min_hindfoot_length max_hindfoot_length
##         <chr>                <dbl>               <int>               <int>
## 1          AH             33.00000                  31                  35
## 2          BA             13.00000                   6                  16
## 3          DM             35.98235                  16                  50
## 4          DO             35.60755                  26                  64
## 5          DS             49.94887                  39                  58
## 6          NL             32.29423                  21                  70
## 7          OL             20.53261                  12                  39
## 8          OT             20.26741                  13                  50
## 9          OX             19.12500                  13                  21
## 10         PB             26.11592                   2                  47
## # ... with 15 more rows
```

```
## Source: local data frame [27 x 4]
## Groups: year [26]
## 
##     year     genus     species weight
##    <int>     <chr>       <chr>  <int>
## 1   1977 Dipodomys spectabilis    149
## 2   1978   Neotoma    albigula    232
## 3   1978   Neotoma    albigula    232
## 4   1979   Neotoma    albigula    274
## 5   1980   Neotoma    albigula    243
## 6   1981   Neotoma    albigula    264
## 7   1982   Neotoma    albigula    252
## 8   1983   Neotoma    albigula    256
## 9   1984   Neotoma    albigula    259
## 10  1985   Neotoma    albigula    225
## # ... with 17 more rows
```

```
## # A tibble: 3 × 2
##     sex     n
##   <chr> <int>
## 1     F 15690
## 2     M 17348
## 3  <NA>  1748
```
--->


See also [`dplyr` cheat sheet](http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)

## Extending tidy data

Note that sometimes, data is not tidy for good reasons - either for
performance reasons, or because of other valid conventions. 


Let's look at a well-known example:


```r
library("Biobase")
data(sample.ExpressionSet)
exprs(sample.ExpressionSet)[1:5, 1:5]
```

```
##                        A         B        C        D        E
## AFFX-MurIL2_at  192.7420  85.75330 176.7570 135.5750 64.49390
## AFFX-MurIL10_at  97.1370 126.19600  77.9216  93.3713 24.39860
## AFFX-MurIL4_at   45.8192   8.83135  33.0632  28.7072  5.94492
## AFFX-MurFAS_at   22.5445   3.60093  14.6883  12.3397 36.86630
## AFFX-BioB-5_at   96.7875  30.43800  46.1271  70.9319 56.17440
```

Typically, slots in the Bioconductor are not (necessarily) tidy. We
could make it tidy with:


```r
library("reshape2")
exprs(sample.ExpressionSet)[1:5, 1:3] %>% melt()
```

```
##               Var1 Var2     value
## 1   AFFX-MurIL2_at    A 192.74200
## 2  AFFX-MurIL10_at    A  97.13700
## 3   AFFX-MurIL4_at    A  45.81920
## 4   AFFX-MurFAS_at    A  22.54450
## 5   AFFX-BioB-5_at    A  96.78750
## 6   AFFX-MurIL2_at    B  85.75330
## 7  AFFX-MurIL10_at    B 126.19600
## 8   AFFX-MurIL4_at    B   8.83135
## 9   AFFX-MurFAS_at    B   3.60093
## 10  AFFX-BioB-5_at    B  30.43800
## 11  AFFX-MurIL2_at    C 176.75700
## 12 AFFX-MurIL10_at    C  77.92160
## 13  AFFX-MurIL4_at    C  33.06320
## 14  AFFX-MurFAS_at    C  14.68830
## 15  AFFX-BioB-5_at    C  46.12710
```

But that wouldn't be helpful for us, working in omics, and could even
counterproductive - see for example the
[Non-tidy data](http://simplystatistics.org/2016/02/17/non-tidy-data/)
blog post.

But we can easily generalise the tidyverse concept:

From tidy data to adequately structured classes:

* Tidy data is data formatted as a table that guarantees that we know
  where to find variables (along columns), observations (along rows),
  and that each cell contains only one value.
  
* For S4 classes, that we use in Bioconductor to store complex data
  that do not fit in rectangular tables, we know where to find every
  bit of information too.

Tools that preserve data tidyness and endomorphism

* The `dplyr` tools take tidy data as input and guarantee to return
  tidy data. 
  
* We should write operations that preserve the classes of the
  data. Ideally, also define simple vocabulary that preserves the rich
  Bioconductor semantic in a consistent paradigm.


```r
library("MSnbase")
library("msdata")
library("magrittr")

fl <- proteomics(full.names = TRUE)[2]
rw <- readMSData2(fl)
rw2 <- rw %>%
    filterMsLevel(3L) %>%
    filterMz(c(126, 132))
```

## References

* [R for Data Science](http://r4ds.had.co.nz/), Hadley Wickham and
  Garrett Grolemund.
* Hadley Wickham, *Tidy Data*, Vol. 59, Issue 10, Sep 2014, Journal of
  Statistical
  Software. [http://www.jstatsoft.org/v59/i10](http://www.jstatsoft.org/v59/i10).
* The tidyverse http://tidyverse.org/, Hadley Wickham.
* [Non-tidy data](http://simplystatistics.org/2016/02/17/non-tidy-data/), Jeff Leek
* The
  [`dplyr` cheat sheet](http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
