context("Testing supporting_function")
library(MetaDICT)
library(testthat)

test_that("Function pcoa_plot_discrete works", {
  data("exampleData")
  O = exampleData$O
  meta = exampleData$meta
  batchid = meta$batch
  p = pcoa_plot_discrete(O,batchid,"Batch")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(p, "ggplot")
    expect_true(length(p$layers) >= 1)
  } else {
    skip("ggplot2 not available")
  }
})


test_that("Function pcoa_plot_continuous works", {
  data("exampleData")
  O = exampleData$O
  Y = runif(ncol(O))
  p = pcoa_plot_continuous(O,Y,"Y")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(p, "ggplot")
    expect_true(length(p$layers) >= 1)
  } else {
    skip("ggplot2 not available")
  }
})

test_that("Function community_detection works", {
  data("exampleData")
  O = exampleData$O
  meta = exampleData$meta
  dist_mat = exampleData$dist_mat
  metadict_res = MetaDICT(O, meta, distance_matrix = dist_mat)
  X = metadict_res$count
  D = metadict_res$D
  D_filter = D[,1:20]
  taxa_c = community_detection(D_filter, max_k = 5)
  expect_type(taxa_c, "list")
})
