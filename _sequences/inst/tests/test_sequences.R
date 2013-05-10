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
  xx <- readFasta(f)
  expect_true(all.equal(xx, dnaseq))
})
