# Task: make your package  "mycamweather"

For this version you can use the .rda file which is a modified version
of the weather-raw.csv file.  See how the .rda file was created below.

You should be able to test the weather file using:

	load("camweatherraw.rda")
	head(camweatherraw)


# Create the skeleton from the camtemp definition.

	package.skeleton("mycamweather", list=c('camtemp'))


# Update the .Rd files

e.g. use the .Rd files I provided.

# Move .rda file into `camweather/data` folder

	mkdir data

	mv ...



# Build, test, install

	R CMD build mycamweather
	R CMD check mycamweather_1.0.tar.gz

	R CMD INSTALL mycamweather_1.0.tar.gz 


# Try it out:

	require(mycamweather)
	data(camweatherraw)
	tail(camweatherraw)
	camtemp("2013-04-01")
	camtemp("2013-08-30")
	camtemp("2013-12-25")

# Add a vignette to the package

Name the vignette `mycamweather-intro.Rmd`.  The vignette should look
something like [this example](
http://htmlpreview.github.com/?https://github.com/lgatto/rbc/blob/master/rr/camweather-files/mycamweather-intro.html)
and be available when you type:

	require(mycamweather)
	vignette("mycamweather-intro")


# Add a test suite to check for leapyears and test dates are valid

e.g. when calling `camtemp(x)`, is x a valid date?  Write a function
`validdate(x)` and add it to the camtemp() function so that the function
won't attempt to print temperature for invalid dates.  (Hint: an
invalid date could be regarded as any date not in the database.)

# Further ideas to build into package

* Consider rebuilding the database when new data is available on
  cl.cam, or default to reading latest datafile from the web if you
  are online.

* R does a good job on handling dates and times.  See
  [R news 4/1 Help desk](http://www.r-project.org/doc/Rnews/Rnews_2004-1.pdf) e.g.

	
		as.Date("2009_03_20", "%Y_%m_%d") + 15

		t1 = Sys.time()
		t2 = Sys.time()

		t2 - t1

* Ask interesting questions with the data, e.g. how "normal" is today
  compared to the seasonal average?  What is the hottest day of a
  year?

* How well does hours of sunshine (or temperature) correlate with
  number of daylight hours in Cambrige?  Use a simple
  [algorithm](http://quantitative-ecology.blogspot.co.uk/2007/10/approximate-sunrise-and-sunset-times.html)
  to estimate sunrise and sunset times.  (Find lat/long of Computer
  lab.)


  

# SJE Notes

## Creating the Data file

Notes for SJE:

The Rda file is 1Mb rather than 16Mb raw.

First, I added a header myself to the .csv.gz file, then re-saved it:

	camweatherraw <- read.csv('camweather/data/camweatherraw.csv.gz', header=TRUE)
	save(camweatherraw, file='camweatherraw.rda', compress="xz")

## Updating the vignette output

Copy the html file across from this location:

	system.file("doc/mycamweather-intro.html", package="mycamweather")
