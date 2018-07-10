context("Defining a genetic system")
library(peas)

test_that("newGenopheno creates genopheno object", {
  expect_equal(class(newGenopheno()), "genopheno")
  expect_equal(is.list(newGenopheno()), TRUE)
  expect_equal(length(newGenopheno()), 2)
  expect_equal(is.list(newGenopheno()[["geno"]]), TRUE)
  expect_equal(is.null(newGenopheno()[["pheno"]]), TRUE)
})
