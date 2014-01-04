## Create a set of tests to evaluate the leap year function.
## 2014-01-02
##
##
## require(testthat)
## test_file('leapyear.R')
##
## You will need to fix the function to get perfect results!
## There are three sets of tests to work through.

## Simple version.
isleapyear <- function(year) {
  ## Return TRUE if YEAR (an integer) is a leap year.
  ## This version is incomplete.
  (year %% 4) == 0
}



context("Current years")
test_that("Simple tests for recent years:", {
  expect_false( isleapyear(2013))
  expect_false( isleapyear(2014))
  expect_false( isleapyear(2015))
  expect_true( isleapyear(2016))
})

## Problem 1. Handle the cases of centuries correctly.
##
context("Centuries")
test_that("Centuries are handled correctly", {
  expect_false( isleapyear(2100))
  expect_true ( isleapyear(2000))
  expect_false( isleapyear(1900))
  expect_true ( isleapyear(1600))
})


## Problem 2. Handle floating point numbers.
## If a number is given, but is not an integer, we round to nearest
## integer using round(x,0), and give a warning.

context("Handling floating point numbers")
test_that("Floating point years", {
  expect_warning( x <- isleapyear(2016.4))
  expect_true   ( x )
})
# Ex: why have we written the previous two tests using the variable x?
# Hint: try and rewrite the tests without the variable.

## Problem 3.
## Some strings can be converted to numbers.  See if we can cope with those
## too.
## Hints: is.character(), as.numeric()

context("Strings can be numbers")
test_that("Strings as years", {
  expect_error( isleapyear("Nonsense"))
  expect_warning( x <- isleapyear("2012.4"))
  #Add following test only when you are working on it, else the
  #output from testing can be quite noisy.
  #expect_true( x )
})
