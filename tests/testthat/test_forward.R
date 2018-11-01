context("Forward")

test_that("forward_sample computes a valid probability", {
  p <- 50
  K <- 5

  p_init <- rep(1 / K, K)
  p_trans <- array(runif((p - 1) * K * K), c(p - 1, K, K))
  for (j in 1:(p - 1)) {
    p_trans[j, , ] <- p_trans[j, , ] / (matrix(rowSums(p_trans[j, , ]), ncol = 1) %*% rep(1, K))
  }
  p_emit <- array(stats::runif(p * 3 * K), c(p, 3, K))
  for (j in 1:p) {
    p_emit[j, , ] <- p_emit[j, , ] / (matrix(rep(1, 3), ncol = 1) %*% colSums(p_emit[j, , ]))
  }

  x <- (runif(p, min = 0, max = 1) < runif(p, min = 0, max = 1)) +
    (runif(p, min = 0, max = 1) < runif(p, min = 0, max = 1))
  log_prob <- forward_sample(x, p_init, p_trans, p_emit)

  expect_lte(exp(log_prob), 1)
  expect_gte(exp(log_prob), 0)
})

test_that("forward concatenates the sample results", {
  p <- 50
  K <- 5

  p_init <- rep(1 / K, K)
  p_trans <- array(runif((p - 1) * K * K), c(p - 1, K, K))
  for (j in 1:(p - 1)) {
    p_trans[j, , ] <- p_trans[j, , ] / (matrix(rowSums(p_trans[j, , ]), ncol = 1) %*% rep(1, K))
  }
  p_emit <- array(stats::runif(p * 3 * K), c(p, 3, K))
  for (j in 1:p) {
    p_emit[j, , ] <- p_emit[j, , ] / (matrix(rep(1, 3), ncol = 1) %*% colSums(p_emit[j, , ]))
  }

  n_samples <- 100
  X <- matrix((runif(n_samples * p, min = 0, max = 1) < runif(n_samples * p, min = 0, max = 1)) +
    (runif(n_samples * p, min = 0, max = 1) < runif(n_samples * p, min = 0, max = 1)),
  ncol = p, nrow = n_samples
  )

  expect_length(forward(X, p_init, p_trans, p_emit, ncores = 1), n_samples)
  expect_length(forward(X, p_init, p_trans, p_emit, ncores = 3), n_samples)
  expect_equal(
    forward(X, p_init, p_trans, p_emit, ncores = 1),
    forward(X, p_init, p_trans, p_emit, ncores = 3)
  )
})

test_that("the output of cond_prob are valid propensity scores regardless of
          the binarization rules", {
  p <- 50
  K <- 5

  p_init <- rep(1 / K, K)
  p_trans <- array(runif((p - 1) * K * K), c(p - 1, K, K))
  for (j in 1:(p - 1)) {
    p_trans[j, , ] <- p_trans[j, , ] / (matrix(rowSums(p_trans[j, , ]), ncol = 1) %*% rep(1, K))
  }
  p_emit <- array(stats::runif(p * 3 * K), c(p, 3, K))
  for (j in 1:p) {
    p_emit[j, , ] <- p_emit[j, , ] / (matrix(rep(1, 3), ncol = 1) %*% colSums(p_emit[j, , ]))
  }

  hmm <- list()
  hmm[["pInit"]] <- p_init
  hmm[["Q"]] <- p_trans
  hmm[["pEmit"]] <- p_emit

  n_samples <- 100
  X <- matrix((runif(n_samples * p, min = 0, max = 1) < runif(n_samples * p, min = 0, max = 1)) +
    (runif(n_samples * p, min = 0, max = 1) < runif(n_samples * p, min = 0, max = 1)),
  ncol = p, nrow = n_samples, dimnames = list(NULL, paste0("SNP_", 1:p))
  )

  expect_equal(rowSums(cond_prob(X, sample(colnames(X), 1), hmm,
    ncores = 3, binary = TRUE
  )), rep(1, n_samples))
  expect_equal(rowSums(cond_prob(X, sample(colnames(X), 1), hmm,
    ncores = 3, binary = FALSE
  )), rep(1, n_samples))
})
