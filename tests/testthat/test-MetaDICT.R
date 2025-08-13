context("Testing MetaDICT")
library(MetaDICT)
library(testthat)

test_that("MetaDICT provides corrected count tables", {
  data(exampleData)
  O = exampleData$O
  meta = exampleData$meta
  dist_mat = exampleData$dist_mat
  metadict_res = MetaDICT(O,meta,distance_matrix = dist_mat,
                         customize_parameter = TRUE, alpha = 0.01, beta = 0.01)
  X = metadict_res$count
  expect_equal(ncol(X),ncol(O))
})
