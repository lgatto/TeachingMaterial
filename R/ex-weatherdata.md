# Weather data

Using the `weatherdata()` function from the `camweather` package,
download a weather data frame of your choice.


```r
library("camweather")
camweather()
```

```
## Cambridge Weather Data
## Package: camweather 0.1.2
## Credit: Digital Technology Group
##         Computer Laboratory
##         University of Cambridge
## Source: http://www.cl.cam.ac.uk/research/dtg/weather/
## - 6662 weather files
## - from: 1995-06-30
## - to: 2014-01-03
```

```r
x <- weatherdata("2013-03-13")
```


- What weather data is available? See `?weatherdata` and inspect the
  data frame's column names.


```r
names(x)
```

```
##  [1] "Time"           "Temp [degC]"    "Humid [%]"      "DewPt [degC]"  
##  [5] "Press [mBar]"   "WindSp [knots]" "WindDr"         "Sun [hours]"   
##  [9] "Rain [mm]"      "Start"          "MxWSpd [knots]" "Day"
```

```r
head(x)
```

```
##                  Time Temp [degC] Humid [%] DewPt [degC] Press [mBar]
## 1 2013-03-13 00:00:00         1.2        81         -1.7         1007
## 2 2013-03-13 00:30:00         0.8        84         -1.6         1007
## 3 2013-03-13 01:00:00         0.4        85         -1.8         1007
## 4 2013-03-13 01:30:00         0.0        86         -2.1         1007
## 5 2013-03-13 02:00:00        -0.4        88         -2.1         1007
## 6 2013-03-13 02:30:00         0.4        87         -1.5         1007
##   WindSp [knots] WindDr Sun [hours] Rain [mm] Start MxWSpd [knots]
## 1            4.0     NW           0         0 00:00             14
## 2            5.6     NW           0         0 00:00             14
## 3            3.5     NW           0         0 00:00             12
## 4            3.1     NW           0         0 00:00              8
## 5            2.3     NW           0         0 00:00              6
## 6            3.5     NW           0         0 00:00              8
##          Day
## 1 2013-03-13
## 2 2013-03-13
## 3 2013-03-13
## 4 2013-03-13
## 5 2013-03-13
## 6 2013-03-13
```



```r
`?`(weatherdata)
```


- What were the highest and lowest temperatures on that day? Hint: see
  `min`, `max` and/or `range` functions.


```r
x <- nounits(x)
range(x$Temp)
```

```
## [1] -1.6  4.8
```

```r
min(x$Temp)
```

```
## [1] -1.6
```

```r
max(x$Temp)
```

```
## [1] 4.8
```

```r
summary(x$Temp)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   -1.60    0.00    0.40    1.07    2.40    4.80
```


- The `Sun` and `Rain` values are cumulative from `Start`. What is the
  average rainfall per hours for that day? Hint: see `diff` for to
  calculate differences between successive values and `mean`.


```r
mean(diff(x$Rain))
```

```
## [1] 0.07458
```

```r
mean(diff(x$Sun))
```

```
## [1] 0.1027
```


- In what direction has the wind blown most on that day? Hint:
  `table`.


```r
table(x$WindDr)
```

```
## 
##  N NW  W 
##  2 42  5
```

