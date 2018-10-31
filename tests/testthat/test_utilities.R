context("Utilities")

test_that("sample_SNP has all of the required fields in the required length", {
  nX <- 5
  nY <- 5
  nZ12 <- 5
  clusters <- rep(1:25, each = 3)
  names(clusters) <- paste0("SNP_", 1:length(clusters))
  thresh_MAF <- 0.2
  MAF <- runif(length(clusters), min = thresh_MAF, max = .5)
  window_size <- 2

  overlap_marg <- 0
  overlap_inter <- 0
  sampled <- sample_SNP(
    nX, nY, nZ12, clusters,
    MAF, thresh_MAF, window_size, overlap_marg, overlap_inter
  )

  expect_setequal(
    names(sampled),
    c("target", "syner", "marginal", "inter1", "inter2")
  )
  expect_equal(length(sampled$target), 1)
  expect_equal(length(sampled$syner), nX)
  expect_equal(length(sampled$marginal), nY)
  expect_equal(length(sampled$inter1), nZ12)
  expect_equal(length(sampled$inter2), nZ12)
})

test_that("sampled_SNP correctly samples in overlapping scenarios", {
  nX <- 5
  nY <- 5
  nZ12 <- 5
  clusters <- rep(1:25, each = 3)
  names(clusters) <- paste0("SNP_", 1:length(clusters))
  thresh_MAF <- 0.2
  MAF <- runif(length(clusters), min = thresh_MAF, max = .5)
  window_size <- 2

  overlap_marg <- 2
  overlap_inter <- 2
  sampled <- sample_SNP(
    nX, nY, nZ12, clusters,
    MAF, thresh_MAF, window_size, overlap_marg, overlap_inter
  )
  expect_equal(length(intersect(sampled$syner, sampled$marginal)), overlap_marg)
  expect_equal(length(intersect(sampled$syner, sampled$inter1)), overlap_inter)
  expect_equal(length(intersect(sampled$marginal, sampled$inter1)), 0)
})

test_that("gen_model output has the right structure", {
  nX <- 5
  nY <- 5
  nZ12 <- 5

  model_coef <- gen_model(nX, nY, nZ12, mean = rep(1, 4), sd = rep(1, 4))
  expect_setequal(names(model_coef), c("syner", "marg", "inter"))
  expect_setequal(names(model_coef$syner), c("A0", "A1"))
})

test_that("gen_model is not affected by the type of the arguments mean and sd", {
  nX <- 5
  nY <- 5
  nZ12 <- 5

  means <- rnorm(4)
  sds <- runif(4, min = .5, max = 1.5)
  set.seed(8487)
  input_scalar <- gen_model(nX, nY, nZ12, mean = means, sd = sds)

  means <- list(rep(means[1], nX), rep(means[2], nX), rep(means[3], nY), rep(means[4], nZ12))
  sds <- list(rep(sds[1], nX), rep(sds[2], nX), rep(sds[3], nY), rep(sds[4], nZ12))
  set.seed(8487)
  input_vector <- gen_model(nX, nY, nZ12, mean = means, sd = sds)

  expect_identical(input_scalar, input_vector)

})

test_that("sim_phenotype is compatible with the whole pipeline", {
  nX <- 5
  nY <- 5
  nZ12 <- 5
  clusters <- rep(1:25, each = 3)
  names(clusters) <- paste0("SNP_", 1:length(clusters))
  thresh_MAF <- 0.2
  MAF <- runif(length(clusters), min = thresh_MAF, max = .5)
  window_size <- 2
  overlap_marg <- 0
  overlap_inter <- 0

  n_samples <- 200
  p_thresh <- 0.4
  X <- matrix((runif(n_samples * length(clusters), min = 0, max = 1) < p_thresh) +
                (runif(n_samples * length(clusters), min = 0, max = 1) < p_thresh),
              ncol = length(clusters), nrow = n_samples
  )
  colnames(X) <- names(clusters)

  causal <- sample_SNP(
    nX, nY, nZ12, clusters,
    MAF, thresh_MAF, window_size, overlap_marg, overlap_inter
  )
  model <- gen_model(nX, nY, nZ12, mean = rnorm(4), sd = runif(4, min = .5, max = 1.5))
  expect_type(sim_phenotype(X, causal, model, intercept = TRUE), "logical")

})

test_that("merge_cluster runs correctly for k as integer ", {
  center <- 3
  clusters <- rep(1:5, each = 2)
  k <- 1
  merged <- clusters
  merged[merged %in% c(3, 4)] <- 2
  merged[merged == 5] <- 3
  expect_equal(merge_cluster(clusters, center, k), merged)

  center <- 5
  merged <- clusters
  merged[merged == 5] <- 4
  expect_equal(merge_cluster(clusters, center, k), merged)
})

test_that("merge_cluster runs correctly for k as vector", {
  center <- 3
  clusters <- rep(1:5, each = 2)
  k <- c(2, 4, 5)
  merged <- clusters
  merged[merged %in% 3:5] <- 2
  expect_equal(merge_cluster(clusters, center, k), merged)
})
