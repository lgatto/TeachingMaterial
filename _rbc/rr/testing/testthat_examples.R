## Here are some example expectations.
##
## require(testthat)
##
## Many tests have a long form and a short-form.
## The long form is more readable (but longer to type!)

## What is the difference between
## source('testthat_examples.R')
## test_file('testthat_examples.R')

expect_that(10, equals(10))
expect_equal(10, 10)
expect_equal(10, 10)

expect_false(1+2==3)


## Bundle relevant tests together for one function into a context.
context("Strings")

test_that("Check string matching", {
  string <- "Testing is fun!"
  expect_match(string, "Testing")
  expect_match(string, "testing") # case-sensitive
  expect_match(string, "T.*ing")
})

