context("Testing MetaDICT_transfer")
library(MetaDICT)
library(testthat)

test_that("Transfer learning is successfully applied", {
  data("exampleData")
  O = exampleData$O
  meta = exampleData$meta
  dist_mat = exampleData$dist_mat
  metadict_res = MetaDICT(O, meta, distance_matrix = dist_mat)
  X = metadict_res$count
  D = metadict_res$D
  R_list = metadict_res$R
  w = metadict_res$w

  data("exampleData_transfer")
  new_data = exampleData_transfer$new_data
  new_meta = exampleData_transfer$new_meta
  new_data_res = metadict_add_new_data(new_data, new_meta, metadict_res)
  new_count = new_data_res$count
  expect_equal(ncol(new_count),ncol(new_data))

})
