context("testing 'sequences' package'")

test_that("dnaseq validity", {
  data(dnaseq)
  expect_true(validObject(dnaseq))
})

test_that("readFasta", {
  ## loading _valid_ dnaseq
  data(dnaseq)
  ## reading fasta sequence
  f <- dir(system.file("extdata",package="sequences"),pattern="fasta",full.names=TRUE)
  xx <- readFasta(f[1])
  expect_true(all.equal(xx, dnaseq))
})

test_that("ccpp code", {
  gccountr <-
    function(x) tabulate(factor(strsplit(x, "")[[1]]))
  x <- "AACGACTACAGCATACTAC"
  expect_true(identical(gccount(x), gccountr(x)))
  expect_true(identical(gccount2(x), gccountr(x)))
})
