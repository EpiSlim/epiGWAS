context("BOOST")

test_that("BOOST runs correctly in multithreaded and
          singlethreaded enivornment", {
  n_samples <- 500
  n_snps <- 200
  p_thresh <- 0.4
  A <- (runif(n_samples, min = 0, max = 1) < p_thresh) + (runif(n_samples,
    min = 0, max = 1
  ) < p_thresh)
  X <- matrix((runif(n_samples * n_snps, min = 0, max = 1) < p_thresh) +
    (runif(n_samples * n_snps, min = 0, max = 1) < p_thresh),
  ncol = n_snps,
  nrow = n_samples
  )
  Y <- (runif(n_samples, min = 0, max = 1) < 0.5)
  expect_type(BOOST(A, X, Y, ncores = 2), "double")
  expect_type(BOOST(A, X, Y, ncores = 1), "double")
  expect_equal(BOOST(A, X, Y, ncores = 2), BOOST(A, X, Y, ncores = 1))
})

test_that("BOOST computed scores are correct", {
  n_samples <- 500
  n_snps <- 200
  p_thresh <- 0.4
  A <- (runif(n_samples, min = 0, max = 1) < p_thresh) + (runif(n_samples,
    min = 0, max = 1
  ) < p_thresh)
  X <- (runif(n_samples, min = 0, max = 1) < p_thresh) + (runif(n_samples,
    min = 0, max = 1
  ) < p_thresh)
  Y <- (runif(n_samples, min = 0, max = 1) < 0.5)
  data_AX <- data.frame(x1 = factor(A), x2 = factor(X), y = Y)
  fit01 <- glm(y ~ x1 + x2, family = "binomial", data = data_AX)
  fit2 <- glm(y ~ x1 + x2 + x1 * x2, family = "binomial", data = data_AX)

  expect_equal(BOOST(A, matrix(X, ncol = 1), Y), fit01$dev - fit2$dev)
})
4
4
