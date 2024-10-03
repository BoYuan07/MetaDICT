context("Testing MetaDICT")
library(MetaDICT)
library(testthat)

test_that("MetaDICT provides corrected count tables", {
  data(exampleData)
  alpha = 0.01
  beta = 1
  gamma = 0.1
  metadict_res = metadict(O_list,dist,meta.list = meta.list,alpha = alpha,beta = beta,gamma = gamma)
  X = metadict_res$X
  O = do.call(cbind,O_list)
  expect_equal(ncol(X),ncol(O))
})
