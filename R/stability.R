#'
#' @param n original sample size
#' @param size subsample size
#' @param n_subsample total number of subsamples
#'
#' @return
#'
subsample <- function(n, size = n %/% 2, n_subsample) {
  idx <- array(NULL, dim = c(size, n_subsample))
  for (i in 1:n_subsample) {
    idx[, i] <- sample.int(n = n, size = size, replace = FALSE)
  }

  return(idx)
}

#'
#'
#' remark about keeping to default
#'
#' @param X_desc
#' @param Y
#' @param family
#' @param n_subsample
#' @param n_lambda
#' @param lambda_min_ratio
#' @param eps regularization parameter.
#' @param short . If \code{TRUE}, the areas under the curve are computed
#' on the first half of the stability path. +
#' @param verbose
#'
#' @return
#'
#' @export
stabilityBIG <- function(X, Y, family, n_subsample = 20, n_lambda = 100,
                         lambda_min_ratio = 0.01, eps = 1e-5, short = TRUE,
                         ncores = 4, dir = "tmp", prefix = "subX") {
  stopifnot(family %in% c("gaussian", "binomial"))
  stopifnot(bigmemory::is.big.matrix(X))
  stopifnot(dir.exists(dir))

  requireNamespace("bigpca", quietly = FALSE)
  requireNamespace("bigmemory", quietly = FALSE)
  requireNamespace("biglasso", quietly = FALSE)

  idx <- subsample(length(Y), size = (length(Y) %/% 2), n_subsample = n_subsample)

  full_fit <- biglasso::biglasso(
    X = X, y = Y,
    penalty = "enet", family = family,
    nlambda = n_lambda, ncores = ncores,
    lambda.min = lambda_min_ratio, alpha = 1 - eps, warn = FALSE
  )

  if (length(full_fit$lambda) < n_lambda) {
    full_fit <- biglasso::biglasso(
      X = X, y = Y,
      penalty = "enet", family = family,
      nlambda = n_lambda, ncores = ncores,
      lambda.min = min(full_fit$lambda) / (0.99 * max(full_fit$lambda)),
      alpha = 1 - eps, warn = FALSE
    )
  }

  length_lambda <- length(full_fit$lambda)
  stab <- array(0, dim = c(length_lambda, dim(X)[2]))

  for (i in 1:n_subsample) {
    message("----- Running stability selection for subsample ", i, " -----")
    sub_X_desc <- bigpca::big.select(
      X,
      select.rows = idx[, i], select.cols = 1:ncol(X),
      delete.existing = TRUE, deepC = TRUE,
      dir = dir, pref = prefix
    )
    sub_X <- bigpca::get.big.matrix(sub_X_desc)
    sub_X_shared <- bigmemory::deepcopy(sub_X, shared = FALSE, type = "double")
    partial_fit <- biglasso::biglasso(
      X = sub_X_shared, y = Y[idx[, i]],
      penalty = "enet", family = family,
      ncores = ncores, lambda = full_fit$lambda, alpha = 1 - eps, warn = FALSE
    )
    if (length(partial_fit$lambda) < length_lambda) length_lambda <- length(partial_fit$lambda)
    partial_coef <- as.matrix(coef(partial_fit, s = full_fit$lambda)[2:(dim(X)[2] + 1), 1:length_lambda])
    stab <- stab[1:length_lambda, ]
    stab <- stab + t(partial_coef != 0)
  }

  unlink(sub_X_desc, force = TRUE)
  rm(sub_X, sub_X_shared)
  gc()

  stab <- stab / n_subsample

  aucs <- sapply(1:ncol(stab), function(i) {
    DescTools::AUC(
      x = 1:((1 - short) * length_lambda + short * (length_lambda %/% 2)),
      y = stab[1:((1 - short) * length_lambda + short * (length_lambda %/% 2)), i]
    )
  })

  return(aucs)
}

#' area under the stability path
#'
#' remark about keeping to default
#'
#' @param X
#' @param Y response vector
#' @param weights
#' @param family
#' @param n_subsample
#' @param n_lambda
#' @param short
#' @param lambda_min_ratio
#' @param eps
#'
#' @return
#'
#' @export
stabilityGLM <- function(X, Y, weights, family, n_subsample = 20, n_lambda = 100,
                         short = TRUE, lambda_min_ratio = 0.01, eps = 1e-5) {
  stopifnot(family %in% c("gaussian", "binomial"))
  stopifnot(all(weights >= 0))

  idx <- subsample(length(Y), size = (length(Y) %/% 2), n_subsample = n_subsample)

  full_fit <- glmnet::glmnet(
    x = X, y = Y, weights = weights,
    family = family, nlambda = n_lambda,
    lambda.min.ratio = 0.01, alpha = 1 - eps
  )

  if (length(full_fit$lambda) < n_lambda) {
    full_fit <- glmnet::glmnet(
      x = X, y = Y, weights = weights,
      family = family, nlambda = n_lambda, alpha = 1 - eps,
      lambda.min.ratio = min(full_fit$lambda) / (0.99 * max(full_fit$lambda))
    )
  }

  stab <- array(0, dim = c(n_lambda, ncol(X)))

  for (i in 1:n_subsample) {
    partial_fit <- glmnet::glmnet(
      x = X[idx[, i], ], y = Y[idx[, i]], weights = weights[idx[, i]],
      family = family, lambda = full_fit$lambda, alpha = 1 - eps
    )

    partial_coef <- glmnet::coef.glmnet(partial_fit, s = full_fit$lambda)[-1, ]
    stab <- stab + (t(as.matrix(partial_coef)) != 0)
  }

  stab <- stab / n_subsample

  aucs <- sapply(1:ncol(X), function(i) {
    DescTools::AUC(
      x = 1:((1 - short) * n_lambda + short * (n_lambda %/% 2)),
      y = stab[1:((1 - short) * n_lambda + short * (n_lambda %/% 2)), i]
    )
  })

  return(aucs)
}
