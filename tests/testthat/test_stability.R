context("Stability")

test_that("subsampling is performed without replacement", {
  n <- 50
  n_subsample <- 10
  sampled <- subsample(n, size = n %/% 2, n_subsample)
  expect_equal(dim(sampled), c(n %/% 2, n_subsample))
  expect_true(all(duplicated.array(sampled, MARGIN = 2) == FALSE))

})
