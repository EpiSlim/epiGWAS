#' Creates multiple subsamples without replacement
#'
#' The subsampling is iteratively performed in order to generate
#' multiple subsamples of a predetermined size.
#'
#' @param n original sample size
#' @param size subsample size
#' @param n_subsample total number of subsamples
#'
#' @return a matrix of indices with \code{size} rows
#' and \code{n_subsample} columns.
#'
#' @examples
#' n <- 50 # Total number of samples
#' n_subsample <- 10 # Number of subsamples
#'
#' sub_matrix <- subsample(n = n, n_subsample = n_subsample)
#' @export
subsample <- function(n, size = n %/% 2, n_subsample) {
  idx <- array(dim = c(size, n_subsample))
  for (i in seq_len(n_subsample)) {
    idx[, i] <- sample.int(n = n, size = size, replace = FALSE)
  }

  return(idx)
}

#' Computes the area under the stability path for all covariates
#'
#' To perform model selection, this function scores all covariates in \code{X}
#' according to the area under their stability selection paths. Our model selection
#' procedure starts by dynamically defining a grid for
#' the elastic net penalization parameter \eqn{\lambda}{lambda}. To define the
#' grid, we solve the full-dataset elastic net. This yields
#' \code{n_lambda} log-scaled values between \eqn{\lambda_{max}}{lambda_max}
#' and \eqn{\lambda_{min}}{lambda_min}. \eqn{\lambda_{max}}{lambda_max} is
#' the maximum value for which the elastic net support is not empty. On the other hand,
#' \eqn{\lambda_{min}}{lambda_min} can be derived through
#' \code{lambda_min_ratio}, which is the ratio of \eqn{\lambda_{min}}{lambda_min}
#' to \eqn{\lambda_{max}}{lambda_max}. The next step is identical to the original
#' stability selection procedure. For each value of \eqn{\lambda}{lambda}, we
#' solve \code{n_subsample} times the same elastic net, though for a different subsample.
#' The subsample is a random selection of half of the samples of the original dataset.
#' The empirical frequency of each covariate entering the support is then the number of
#' times the covariate is selected in the support as a fraction of \code{n_subsample}.
#' We obtain the stability path by associating each value of \eqn{\lambda}{lambda}
#' with the corresponding empirical frequency. The final scores are the areas under the
#' stability path curves. This is a key difference with the original stability
#' selection procedure where the final score is the maximum empirical frequency.
#' On simulations, our scoring technique outperformed maximum empirical
#' frequencies.
#'
#' @family support estimation functions
#'
#' @param X input design matrix
#' @param Y response vector
#' @param weights nonnegative sample weights
#' @param family response type. Either 'gaussian' or 'binomial'
#' @param n_subsample number of subsamples for stability selection
#' @param n_lambda total number of lambda values
#' @param short whether to compute the aucs only on the first half
#' of the stability path. We observed better performance for
#' thresholded paths
#' @param lambda_min_ratio ratio of \eqn{\lambda_{min}}{lambda_min} to
#' \eqn{\lambda_{max}}{lambda_max} (see description for a thorough explanation)
#' @param eps elastic net mixing parameter.
#'
#' @return a vector containing the areas under the stability path
#' curves
#'
#' @details For a fixed \eqn{\lambda}{lambda},
#' the L2 penalization is \eqn{\lambda \times eps}{lambda * eps}, while
#' the L1 penalization is \eqn{\lambda \times (1-eps)}{lambda * (1-eps)}.
#' The goal of the L2 penalization is to ensure the uniqueness of the
#' solution. For that reason, we recommend setting eps << 1.
#'
#' @references Slim, L., Chatelain, C., Azencott, C.-A., & Vert, J.-P.
#' (2018). Novel Methods for Epistasis Detection in Genome-Wide Association
#' Studies. BioRxiv.
#'
#' @references Meinshausen, N., & Bühlmann, P. (2010). Stability
#' selection. Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology), 72(4), 417–473.
#'
#' @references Haury, A. C., Mordelet, F., Vera-Licona, P., & Vert, J. P.
#' (2012). TIGRESS: Trustful Inference of Gene REgulation using Stability
#' Selection. BMC Systems Biology, 6.
#'
#' @examples
#' # ---- Continuous data ----
#' n <- 50
#' p <- 20
#' X <- matrix(rnorm(n * p), ncol = p)
#' Y <- crossprod(t(X), rnorm(p))
#' aucs_cont <- stabilityGLM(X, Y,
#'   family = "gaussian", n_subsample = 1,
#'   short = FALSE
#' )
#'
#' # ---- Binary data ----
#' X <- matrix(rnorm(n * p), ncol = p)
#' Y <- runif(n, min = 0, max = 1) < 1 / (1 + exp(-X[, c(1, 7, 15)] %*% rnorm(3)))
#' weights <- runif(n, min = 0.4, max = 0.8)
#' aucs_binary <- stabilityGLM(X, Y,
#'   weights = weights,
#'   n_lambda = 50, lambda_min_ratio = 0.05, n_subsample = 1
#' )
#' @seealso \code{\link[glmnet]{glmnet-package}}
#'
#' @export
stabilityGLM <- function(X, Y, weights = rep(1, nrow(X)), family = "gaussian",
                         n_subsample = 20, n_lambda = 100, short = TRUE, lambda_min_ratio = 0.01,
                         eps = 1e-05) {
  stopifnot(family %in% c("gaussian", "binomial"))
  stopifnot(all(weights >= 0))

  idx <- subsample(length(Y), size = (length(Y) %/% 2), n_subsample = n_subsample)

  full_fit <- glmnet::glmnet(
    x = X, y = Y, weights = weights, family = family,
    nlambda = n_lambda, lambda.min.ratio = 0.01, alpha = 1 - eps
  )

  if (length(full_fit$lambda) < n_lambda) {
    full_fit <- glmnet::glmnet(
      x = X, y = Y, weights = weights, family = family,
      nlambda = n_lambda, alpha = 1 - eps, lambda.min.ratio = min(full_fit$lambda) / (0.99 *
        max(full_fit$lambda))
    )
  }

  length_lambda <- length(full_fit$lambda)
  stab <- array(0, dim = c(length_lambda, dim(X)[2]))

  for (i in seq_len(n_subsample)) {
    partial_fit <- glmnet::glmnet(
      x = X[idx[, i], ], y = Y[idx[, i]],
      weights = weights[idx[, i]], family = family, lambda = full_fit$lambda,
      alpha = 1 - eps
    )
    if (length(partial_fit$lambda) < length_lambda) {
      length_lambda <- length(partial_fit$lambda)
    }
    partial_coef <- glmnet::coef.glmnet(partial_fit, s = full_fit$lambda)[seq_len(dim(X)[2]) +
      1, seq_len(length_lambda)]
    stab <- stab[seq_len(length_lambda), ]
    stab <- stab + (t(as.matrix(partial_coef)) != 0)
  }

  stab <- stab / n_subsample

  aucs <- vapply(seq_len(ncol(X)), function(i) {
    DescTools::AUC(
      x = seq_len((1 - short) * n_lambda + short * (n_lambda %/% 2)),
      y = stab[
        seq_len((1 - short) * n_lambda + short * (n_lambda %/% 2)),
        i
      ]
    )
  }, double(1))

  return(aucs)
}


#' Computes the area under the stability path for all covariates
#'
#' This function implements the same model selection technique extensively
#' described in \code{\link{stabilityGLM}}. The sole difference is the use
#' of a different elastic net solver. In this function, we make use of
#' \code{\link[biglasso]{biglasso}}. Thanks to its parallel
#' backend, \code{biglasso} scales well to
#' high-dimensional GWAS datasets. However, in our case, because of the use of
#' additional backend files, a slight decrease in runtime is to be expected,
#' compared with \code{\link{stabilityGLM}}.
#'
#' @family support estimation functions
#'
#' @param X design matrix formatted as a
#' \code{\link[bigmemory]{big.matrix}} object
#' @param Y response vector
#' @param family response type. Either 'gaussian' or 'binomial'
#' @param n_subsample number of subsamples for stability selection
#' @param n_lambda total number of lambda values
#' @param lambda_min_ratio the minimum value of the regularization
#' parameter lambda as a fraction of the maximum lambda, the first
#' value for which the elastic net support is not empty.
#' @param eps elastic net mixing parameter (see \code{\link{stabilityGLM}}
#' for more details)
#' @param short whether to compute the aucs only on the first half
#' of the stability path. We observed better performance with
#' thresholded paths
#' @param ncores number of cores for the
#'  \code{\link[biglasso]{biglasso}} solver
#'
#' @return a vector grouping the aucs of all covariates within \code{X}
#'
#' @references Slim, L., Chatelain, C., Azencott, C.-A., & Vert, J.-P.
#' (2018). Novel Methods for Epistasis Detection in Genome-Wide Association
#' Studies. BioRxiv.
#'
#' @references Meinshausen, N., & Bühlmann, P. (2010). Stability
#' selection. Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology), 72(4), 417–473.
#'
#' @references Haury, A. C., Mordelet, F., Vera-Licona, P., & Vert, J. P.
#' (2012). TIGRESS: Trustful Inference of Gene REgulation using Stability
#' Selection. BMC Systems Biology, 6.
#'
#' @seealso \code{\link[biglasso]{biglasso-package}}
#'
#' @examples
#' n <- 100
#' p <- 25
#' X <- bigmemory::as.big.matrix(matrix(runif(n * p), ncol = p))
#' Y <- runif(n, min = 0, max = 1) < 0.5
#' aucBIG <- stabilityBIG(X, Y,
#'   family = "binomial", short = TRUE,
#'   ncores = 1, n_lambda = 200, n_subsample = 1
#' )
#' @export
stabilityBIG <- function(X, Y, family = "gaussian", n_subsample = 20, n_lambda = 100,
                         lambda_min_ratio = 0.01, eps = 1e-05, short = TRUE, ncores = 2) {
  if (requireNamespace("bigmemory", quietly = TRUE) & requireNamespace("biglasso", quietly = TRUE)) {
    stopifnot(family %in% c("gaussian", "binomial"))
    stopifnot(bigmemory::is.big.matrix(X))

    idx <- subsample(length(Y), size = (length(Y) %/% 2), n_subsample = n_subsample)

    full_fit <- biglasso::biglasso(
      X = X, y = Y, penalty = "enet", family = family,
      nlambda = n_lambda, ncores = ncores, lambda.min = lambda_min_ratio,
      alpha = 1 - eps, warn = FALSE
    )

    if (length(full_fit$lambda) < n_lambda) {
      full_fit <- biglasso::biglasso(
        X = X, y = Y, penalty = "enet", family = family,
        nlambda = n_lambda, ncores = ncores, lambda.min = min(full_fit$lambda) / (0.99 *
          max(full_fit$lambda)), alpha = 1 - eps, warn = FALSE
      )
    }

    length_lambda <- length(full_fit$lambda)
    stab <- array(0, dim = c(length_lambda, dim(X)[2]))

    for (i in seq_len(n_subsample)) {
      sub_X_shared <- bigmemory::deepcopy(X, rows = idx[, i], shared = FALSE, type = "double")
      partial_fit <- biglasso::biglasso(
        X = sub_X_shared, y = Y[idx[, i]],
        penalty = "enet", family = family, ncores = ncores, lambda = full_fit$lambda,
        alpha = 1 - eps, warn = FALSE
      )
      if (length(partial_fit$lambda) < length_lambda) {
        length_lambda <- length(partial_fit$lambda)
      }
      partial_coef <- as.matrix(stats::coef(partial_fit, s = full_fit$lambda)[seq_len(dim(X)[2]) +
        1, seq_len(length_lambda)])
      stab <- stab[seq_len(length_lambda), ]
      stab <- stab + t(partial_coef != 0)
    }
  }

  if (dir.exists("/tmp/boost_interprocess") & grepl("^darwin", R.version$os)) {
    unlink("/tmp/boost_interprocess", recursive = TRUE, force = TRUE) # temporary folder created by biglasso
  }

  remove(sub_X_shared)
  gc()

  stab <- stab / n_subsample

  aucs <- vapply(seq_len(ncol(stab)), function(i) {
    DescTools::AUC(x = seq_len((1 - short) * length_lambda + short *
      (length_lambda %/% 2)), y = stab[seq_len((1 - short) * length_lambda +
      short * (length_lambda %/% 2)), i])
  }, double(1))

  return(aucs)
}
