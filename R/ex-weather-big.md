## Exercise

- Write your own `weatherdata` function that takes a date character as
  input, locates the file in the `camweather` package directory, loads
  it and returns an appropriate data structure. See section on reading
  data for help.


```r
getweather <- function(dt) {
    require("camweather")
    weatherdata(dt)
}
```


- Write a function that takes a date character of the for
  `"YYYY-MM-DD"` as input and produces a plot of temperature over
  time. Make sure that it remains possible to fully customise the
  figure as would be with `plot`.


```r
plotweather <- function(dt, ...) {
    require("camweather")
    x <- nounits(weatherdata(dt))
    plot(x$Time, x$Temp, ...)
}
```


- Select all the weather files from June 2013. All file names are
  available with the `weatherfiles()` function. You can use the `grep`
  function to select the relevant file names. Check that you obtain 30
  files.


```r
fls <- weatherfiles()
```

```
## Error: could not find function "weatherfiles"
```

```r
f <- grep("2013_06", fls, value = TRUE)
```

```
## Error: object 'fls' not found
```

```r
length(f)
```

```
## Error: object 'f' not found
```


- Load the 30 data frames into a convenient data structure. Check the
  number of data points that are available for each weather data set.


```r
xx <- lapply(f, weatherdata)
```

```
## Error: object 'weatherdata' not found
```

```r
sapply(xx, nrow)
```

```
## Error: object 'xx' not found
```

```r
table(sapply(xx, nrow))
```

```
## Error: object 'xx' not found
```


- Calculate the average day temperatures for that month.


```r
sapply(xx, function(x) mean(x[, "Temp [degC]"]))
```

```
## Error: object 'xx' not found
```

```r
## or
dd <- do.call(rbind, xx)
```

```
## Error: object 'xx' not found
```

```r
tapply(dd[, 2], dd$Day, mean)
```

```
## Error: object 'dd' not found
```


- Plot the temperature over the full month and the daily
  temperature curves for June 2013.


```r
plot(dd[, 1], dd[, 2], type = "l",
     xlab = "Time", ylab = "Temp",
     main = "June 2013")
```

```
## Error: object 'dd' not found
```



```r
updateday <- function(x)
    as.POSIXct(strftime(x, "%H:%M"), format = "%H:%M")

library("RColorBrewer")
col <- brewer.pal(10, "Set3")
col[2] <- "#555555"
col <- rep(col, each = 3)
lty <- rep(1:3, 30)

trng <- range(lapply(xx, function(x) x[, "Temp [degC]"]))
```

```
## Error: object 'xx' not found
```

```r
plot(updateday(xx[[1]][, 1]),
     xx[[1]][, 2], ylim = trng, type = "l",
     col = col[1], lty = lty[1], lwd = 2,
     xlab = "Time", ylab = "Temp")
```

```
## Error: object 'xx' not found
```

```r

for (i in 2:length(xx))
    lines(updateday(xx[[i]][, 1]), xx[[i]][, 2],
          col = col[i], lty = lty[i], lwd = 2)
```

```
## Error: object 'xx' not found
```

```r

legend("bottomright", legend = 1:30,
       col = col, lty = lty, lwd = 2,
       bty = "n", cex = .8,
       ncol = 5)
```

```
## Error: plot.new has not been called yet
```

