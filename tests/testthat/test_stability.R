context("Stability")

test_that("subsampling is performed without replacement", {
  n <- 50
  n_subsample <- 10
  sampled <- subsample(n, size = n %/% 2, n_subsample)
  expect_equal(dim(sampled), c(n %/% 2, n_subsample))
  expect_true(all(duplicated.array(sampled, MARGIN = 2) == FALSE))
})

test_that("stabilityGLM correctly adapts to the length of the LASSO path", {
  n <- 200
  p <- 50
  X <- matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  Y <- rnorm(n)
  weights <- runif(n)

  set.seed(468)
  aucs_short <- stabilityGLM(X, Y, weights, family = "gaussian", short = TRUE)
  set.seed(468)
  aucs_long <- stabilityGLM(X, Y, weights, family = "gaussian", short = FALSE)

  expect_true(all(aucs_short <= aucs_long))
})

test_that("stabilityGLM shortens the path if the model saturates", {
  n <- 200
  p <- 50
  X <- matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  Y <- rnorm(n)
  weights <- runif(n)

  expect_type(stabilityGLM(X, Y, weights,
                           family = "gaussian",
                           lambda_min_ratio = 1e-5, short = TRUE
  ), "double")
})

test_that("stabilityBIG correctly adapts to the length of the LASSO path", {
  n <- 200
  p <- 50
  X <- bigmemory::as.big.matrix(
    matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  )
  Y <- rnorm(n)

  set.seed(468)
  aucs_short <- stabilityBIG(X, Y,family = "gaussian", ncores = 3, short = TRUE)
  set.seed(468)
  aucs_long <- stabilityBIG(X, Y, family = "gaussian", ncores = 3, short = FALSE)

  expect_true(all(aucs_short <= aucs_long))
})

test_that("stabilityBIG shortens the path if the model saturates", {
  n <- 200
  p <- 50
  X <- bigmemory::as.big.matrix(
    matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  )
  Y <- rnorm(n)
  weights <- runif(n)

  expect_type(stabilityBIG(X, Y, family = "gaussian", ncores = 3,
                           lambda_min_ratio = 1e-5, short = TRUE
  ), "double")
})
