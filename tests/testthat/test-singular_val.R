context("Testing plot_singular_values")
library(MetaDICT)
library(testthat)

test_that("Function for singular value plots works", {
  data("exampleData")
  O = exampleData$O
  meta = exampleData$meta
  p = plot_singular_values(O, meta)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(p, "ggplot")
    expect_true(length(p$layers) >= 1)
  } else {
    skip("ggplot2 not available")
  }
})
