context("Methods")

test_that("OWL executes the do.call instruction", {
  n_samples <- 200
  p <- 50
  X <- matrix((runif(n_samples * p, min = 0, max = 1) < runif(n_samples * p, min = 0, max = 1)) +
                (runif(n_samples * p, min = 0, max = 1) < runif(n_samples * p, min = 0, max = 1)),
              ncol = p, nrow = n_samples
  )
  A <- (runif(n_samples, min = 0, max = 1) < 0.5)
  propensity <- runif(n_samples, min = 0.4, max = 0.8)
  Y <- (runif(n_samples, min = 0, max = 1) < 0.5)
  expect_length(OWL(A, X, Y, propensity), p)
})
