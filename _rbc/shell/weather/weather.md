% Analyse the Cambridge weather data.
% Stephen Eglen
% 2013-12-31

<!--pandoc
format: html
s:
mathjax:
number-sections:

format: latex
number-sections:
-->


# Show of hands

Who has never used Unix?

# Prerequsites.

The command line

Redirection

What's a pipe?

cat head, tail

basic editors: gedit, TextEdit



## man pages
e.g. man cut

#  Intro: what is it?

[CL.cam weather data](https://www.cl.cam.ac.uk/research/dtg/weather/)

e.g. [average weather in April 2013](https://www.cl.cam.ac.uk/research/dtg/weather/period-graph.cgi?2013-04)

or a [daily graph](https://www.cl.cam.ac.uk/research/dtg/weather/daily-graph.cgi?2013-04-01)

Can we do something like this ourselves with the data?


## Get the daily files.

We could also work from the raw CSV, but let's use these 'daily' files for
now.


```bash
mkdir -p weather
cd weather
wget https://www.cl.cam.ac.uk/research/dtg/weather/weather.tar.gz
tar zxf weather.tar.gz
```



How many files did I just download?

```bash
cd daily-text
ls
ls | wc
```


# Get the temperature for one day.

Save the temp for one day in a file.  We need to extract one column:


```bash
tail -n +9 2013_04_01  ## Argh, let's keep the text.
tail -n +9 2013_04_01 > tmp
ls -l tmp
head tmp
## TAB delimeted.
cut -f1,2 tmp > 2013_04_01-temp # EX: shorten
```



# Plot the weather

What we need is a little prog that takes our daily temperature file
and makes a nice plot of it.


```r
dat <- read.table("2013_04_01-temp")
plot(dat)
plot(dat, ylim = c(-10, 30))
```


So far, so good, but those times don't look right.  Let's change the
colons to periods.



```bash
cut -f1,2 tmp | sed 's/:/./' > 2013_04_01-temp # EX: shorten
```


# Running R in BATCH.

(Did Laurent cover this?)


Ex: Show how to do R CMD BATCH, and where output goes.  Why is this
useful?
<!-- ```{r prelim, engine='cat', engine.opts=list(file='trig.R')} -->

```r
## Example script file, to generate a little plot.
## Demonstrate trigonometric functions.
## Sept 2007
x <- seq(from=0, to=2*pi, length=100)
y <- sin(x)
z <- cos(2*x)
z                               ## will not appear
print(y[1:10])                        ## should use print()
plot(x, y, type='l')
lines(x, z, type='l', col='red')
```


- At the command line, type `R CMD BATCH trig.R`.  R will start up,
  process your commands and then quit. 
- Output is stored in the file `trig.Rout`
- If there were no errors, the last line of the output file
    shows the time taken to run the script.
- Any output is not shown on the screen but sent to a PDF
  called `Rplots.pdf`

- This is a GREAT way of testing your scripts, since R starts with an
  empty workspace, you will see if you have all the steps needed.
    
- Aim to always leave your scripts in a working state at the end of a
  session, so that a few days later you don't have to remember why it
  wasn't working!


# Calling R from the command line.


```bash
Rscript -e 'sqrt(1:40)'
Rscript -e 'round(runif(10))'
```

```
##  [1] 1.000000 1.414214 1.732051 2.000000 2.236068 2.449490 2.645751 2.828427
##  [9] 3.000000 3.162278 3.316625 3.464102 3.605551 3.741657 3.872983 4.000000
## [17] 4.123106 4.242641 4.358899 4.472136 4.582576 4.690416 4.795832 4.898979
## [25] 5.000000 5.099020 5.196152 5.291503 5.385165 5.477226 5.567764 5.656854
## [33] 5.744563 5.830952 5.916080 6.000000 6.082763 6.164414 6.244998 6.324555
##  [1] 1 1 1 1 1 0 1 1 1 1
```


## Passing arguments to an Rscript.

Ex: Create a new file simple_rnorm.R, with contents:

<!-- ```{r, engine='cat', engine.opts=list(file='simple_rnorm.R')} -->


```r
#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

stopifnot(length(args)==3)
args = as.numeric(args)
n = args[1]; mean = args[2]; sd = args[3]
rnorm(n, mean, sd)
```



```bash
chmod +x simple_rnorm.R
./simple_rnorm.R 5 10 2
./simple_rnorm.R 100 
```


EX: update simple_rnorm.R so that mean and s.d. take sensible defaults
if not given on command line.


## Now let's create our weather pdf.

Save this to `maketemp.R` in the same folder as the data.
<!--   ```{r, engine='cat', engine.opts=list(file='maketemp.R'), echo=TRUE} -->


```r
#!/usr/bin/env Rscript
# maketemp.R (master in weather.Rmd)
args <- commandArgs(TRUE)
stopifnot(length(args)==1)
file <- args[1]
dat <- read.table(file)


## Create an output file name based on the input file name
pdfname <- paste0(file, ".pdf")
pdf(file=pdfname)
plot(dat, ylim=c(-10, 30), xlab='Time (h)', ylab='Temp (C)',
     type='l', col='red')
abline(h=0, col='grey')
title(file)
dev.off()
```



# Spare Exercises

Find the hottest (or coldest) temp recorded in the archive.  Can you
do this with unix calls?

On how many days did the tempearature go above 30 C?  Below -10 C?

What was the sunniest (or cloudiest) day in the archive?

When was their the biggest difference in temperature from one day to
the next?  Plot the average temp per day and then find the largest
difference?  When was there the biggest "spike" (e.g. one day there
was a sudden change in temperature)?

Can you write a general script to plot all the variables in the
database?  Our all graphs on one page?  (hint: `par(mfrow=c(3,2)`)).

Are there differences between the content in the daily files and
the raw.csv file?

# Principles
Separate computation from plotting

keep data in a separate read-only folder

Carefully delineate what must be kept vs what can be regenerated.


## End: how to compile.

```r
require(knitr)
knit("weather.Rmd")
knit2html("weather.Rmd")
browseURL("weather.html")
```

