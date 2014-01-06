# Weather plot

Using a weather data frame as input, generate a plot showing the
hourly (or half-hourly) rainfall for the 3rd Jan 2014.


```r
x <- nounits(weatherdata("2014-01-03"))
```

```
## Error: could not find function "nounits"
```

```r
rain <- x$Rain
```

```
## Error: object 'x' not found
```

```r
plot(x$Time[-1], diff(rain), type = "l", xlab = "Time", ylab = "Rain", main = "2014-01-03", 
    col = "blue")
```

```
## Error: object 'x' not found
```

