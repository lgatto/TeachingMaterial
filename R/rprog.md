R programming
===




# IO

## Reading text spreadsheets

Using `read.table` *et al.*, a text spreadsheet (in `csv`, `tsv` or
similar) can be read in and converted into a `data.frame`. As usual,
text columns as converted into factors unless `stringsAsFactor=FALSE`
or `as.is=TRUE`.


```r
f <- "2014_01_01"
lf <- file.path("../data/daily-text", f)
w <- read.table(lf, header = FALSE, comment.char = "#", sep = "\t")
dim(w)
head(w)
```


We are missing header information. Let's use `readLines` to parse and
extract the header data manually. 


```r
hd <- readLines("../data/daily-text/2014_01_01")
hd <- hd[grep("#", hd)]
hd <- sub("#", "", hd)
hd <- hd[7:8]
hd <- gsub(" ", "", hd)
hd <- strsplit(hd, "\t")
hd <- paste0(hd[[1]], " [", hd[[2]], "]")
hd <- sub(" \\[\\]", "", hd)
names(w) <- hd
```


The format of the first columns `Time` is still unsatisfactory. Below,
we convert it into a time format/date.


```r
class(w$Time)
w$Time <- strptime(paste(f, w$Time), "%Y_%m_%d %H:%M")
class(w$Time)
summary(w)
```


## Basic plotting


```r
par(mfrow = c(2, 2))
plot(w$Time, w[, "Temp [degC]"], type = "b", xlab = "Time", ylab = "Temp")
plot(w$Time, w[, "WindSp [knots]"], type = "b", xlab = "Time", ylab = "Wind speed")
plot(w$Time, w[, "Rain [mm]"], type = "b", xlab = "Time", ylab = "Rain")
plot(w$Time, w[, "Press [mBar]"], type = "b", xlab = "Time", ylab = "Pressure")
```



```r
boxplot(w[, "WindSp [knots]"] ~ factor(w$WindDr))
pairs(w[, c(2, 5, 6, 9)])
```


## More plotting

While using two axes can be very misleading when the scales are
different (as in the example below) and the differences are not
properly accounted for, let's illustrate such an example to learn how
to set different elements of a base plot. 

- Data rescaling


```r
temp0 <- w[, "Temp [degC]"]
temp <- temp0 - min(temp0)  ## min is 0
temp <- temp/max(temp)  ## max is 1

press0 <- w[, "Press [mBar]"]
press <- press0 - min(press0)
press <- press/max(press)
```


- Plot with minimal decoration


```r
par(mar = c(5, 4, 2, 4))
plot(w$Time, temp, type = "l", xlab = "Time", ylab = "Temp [deg C]", yaxt = "n", 
    col = "steelblue")
lines(w$Time, press, col = "red")
```


- Axis, title and legends


```r
axis(2, at = seq(0, 1, length = 11), labels = seq(min(temp0), max(temp0), length = 11))
axis(4, at = seq(0, 1, length = 11), labels = seq(min(press0), max(press0), 
    length = 11))
mtext("Pressure [mBar]", 4, line = 3)
title(f)
legend("top", c("Temperature", "Pressure"), col = c("steelblue", "red"), lty = 1, 
    bty = "n")
```


## Writing text spreadsheets

`write.table` and `writeLines` to write tables and lines to files.

## Saving and loading R data

Once data has been loaded an properly formatted into `R`, the easiest
and fastest way to serialise it is with `save`. Once stored in binary
format, the saved variable can be loaded in the working environment
with `load`. This can be used for one, several or all object in the
workspace (see `save.image`).

To save and restore single object to a file, see also `saveRDS` and
`readRDS`.

### Exercise

Use one of the `save`/`load` or `saveRDS`/`readRDS` above to serialise
the weather `data.frame`, rename the variable in your workspace,
reload the saved object and compare it to the original file. Hint: use
`identical` or `all.equal`.

### Bonus


```r
x <- sqrt(2)
x^2 == 2
all.equal(x^2, 2)
```

See [Why doesn’t R think these numbers are equal?](http://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f).

## See also
- `scan`: Read data into a vector or list from the console or file.
- User input: `menu`


```r
choices <- c("RStudio", "Wordpad", "emacs", "vim", "Notepad++")
mychoice <- menu(choices, graphics = FALSE, title = "Best editor ever")
cat("Best editor ever is ", choices[mychoice], "\n")
```

- Database: `RMySQL`, `RMongoDB`, `ROracle`
- `rhdf5`, `ncdf`, `XML`, `RJSONIO`, `jsonlite`
- web: `RCurl`, `httr`

# Iteration and flow control 

- `for` and `while`
- `if`, `else`
- `*apply`
- `plyr`, `reshape2`
- `parallel`
- `replicate`
- `with`

# Writing function
- `function`
- pass by value (vs by reference)
- `environment`
- scoping

# Documentation

- `roxygen2` package syntax for in-line documentation

# Misc

- string processing: `strsplit`, `sub`, `gsub`, `paste`, `cat`
- more regexp: `stringr`, `tm`
- `message`, `warning`, `error`
- timing, benchmarking
- debugging
