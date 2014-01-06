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
something like [this example](mycamweather-intro.html) and be
available when you type:

	require(mycamweather)
	vignette("mycamweather-intro")


# Add a test suite to check for leapyears and test dates are valid

e.g. when calling `camtemp(x)`, is x a valid date?  Write a function
`validdate(x)` and add it to the camtemp() function so that the function
won't attempt to print temperature for invalid dates.  (Hint: an
invalid date could be regarded as any date not in the database.)

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
