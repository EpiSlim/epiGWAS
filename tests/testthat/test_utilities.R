context("Utilities")

# sample_SNP
# gen_model
# sim_phenotype

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
  expect_setequal(length(sampled$target), 1)
  expect_setequal(length(sampled$syner), nX)
  expect_setequal(length(sampled$marginal), nY)
  expect_setequal(length(sampled$inter1), nZ12)
  expect_setequal(length(sampled$inter2), nZ12)
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
  expect_setequal(length(intersect(sampled$syner, sampled$marginal)), overlap_marg)
  expect_setequal(length(intersect(sampled$syner, sampled$inter1)), overlap_inter)
  expect_setequal(length(intersect(sampled$marginal, sampled$inter1)), 0)
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
