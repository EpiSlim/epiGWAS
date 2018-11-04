context("Methods")

test_that("OWL executes the do.call instruction", {
  n_samples <- 200
  p <- 20
  X <- matrix((runif(n_samples * p, min = 0, max = 1) < runif(n_samples *
    p, min = 0, max = 1)) + (runif(n_samples * p, min = 0, max = 1) <
    runif(n_samples * p, min = 0, max = 1)), ncol = p, nrow = n_samples)
  A <- (runif(n_samples, min = 0, max = 1) < 0.5)
  propensity <- runif(n_samples, min = 0.4, max = 0.8)
  Y <- (runif(n_samples, min = 0, max = 1) < 0.5)
  expect_length(OWL(A, X, Y, propensity), p)
  expect_error(OWL(A, X, Y, propensity, family = "gaussian"))
})

test_that("do.call is instructed for stabilityGLM", {
  stability_mode <- FALSE
  n_samples <- 200
  p <- 20
  X <- matrix((runif(n_samples * p, min = 0, max = 1) < runif(n_samples *
    p, min = 0, max = 1)) + (runif(n_samples * p, min = 0, max = 1) <
    runif(n_samples * p, min = 0, max = 1)), ncol = p, nrow = n_samples)
  A <- (runif(n_samples, min = 0, max = 1) < 0.5)
  propensity <- runif(n_samples, min = 0.4, max = 0.8)
  Y <- (runif(n_samples, min = 0, max = 1) < 0.5)

  expect_length(modified_outcome(A, X, Y, propensity,
    parallel = stability_mode,
    n_subsample = 20, n_lambda = 100
  ), p)
  expect_length(shifted_outcome(A, X, Y, propensity,
    shift = 0.2, parallel = stability_mode,
    n_subsample = 20, n_lambda = 100
  ), p)
  expect_length(normalized_outcome(A, X, Y, propensity,
    parallel = stability_mode,
    n_subsample = 20, n_lambda = 100
  ), p)
  expect_length(robust_outcome(A, X, Y, propensity,
    parallel = stability_mode,
    n_subsample = 20, n_lambda = 100
  ), p)
})

test_that("do.call is instructed for stabilityBIG", {
  stability_mode <- TRUE
  n_samples <- 200
  p <- 20
  X <- matrix((runif(n_samples * p, min = 0, max = 1) < runif(n_samples *
    p, min = 0, max = 1)) + (runif(n_samples * p, min = 0, max = 1) <
    runif(n_samples * p, min = 0, max = 1)), ncol = p, nrow = n_samples)
  A <- (runif(n_samples, min = 0, max = 1) < 0.5)
  propensity <- runif(n_samples, min = 0.4, max = 0.8)
  Y <- (runif(n_samples, min = 0, max = 1) < 0.5)

  expect_length(modified_outcome(A, X, Y, propensity,
    parallel = stability_mode,
    ncores = 2, n_subsample = 20, n_lambda = 100
  ), p)
  expect_length(normalized_outcome(A, X, Y, propensity,
    parallel = stability_mode,
    ncores = 2, n_subsample = 20, n_lambda = 100
  ), p)
  expect_length(shifted_outcome(A, X, Y, propensity,
    shift = 0.2, parallel = stability_mode,
    ncores = 2, n_subsample = 20, n_lambda = 100
  ), p)
  expect_length(robust_outcome(A, X, Y, propensity,
    parallel = stability_mode,
    ncores = 2, n_subsample = 20, n_lambda = 100
  ), p)
})
