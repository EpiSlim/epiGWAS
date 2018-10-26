#' Implements the outcome weighted learning approach
#'
#' briefly explain modified outcome
#'
#' @param A target variant. The variable A must be encoded as (0, 1)
#' or (0, 1, 2)
#' @param X rest of the genotype
#' @param Y phenotype (binary or continuous)
#' @param propensity propensity scores (a vector or a two-column matrix)
#' @param ... additional arguments to \code{stabilityGLM}
#'
#' @details If the outcome \code{Y} is not nonnegative, we translate it
#' to get
#'
#' @references Zhao, Y., Zeng, D., Rush, A. J., & Kosorok, M. R. (2012).
#' Estimating Individualized Treatment Rules Using Outcome Weighted Learning.
#' Journal of the American Statistical Association, 107(499), 1106–1118.
#'
#' @export
OWL <- function(A, X, Y, propensity, ...) {
  owl_args <- list(...)
  if ("family" %in% names(owl_args)) {
    stopifnot(owl_args["family"] == "binomial")
  } else {
    owl_args["family"] <- "binomial"
  }

  stopifnot(all(propensity > 0) & (propensity <= 1))
  if (is.matrix(propensity)) stopifnot(dim(propensity)[2] == 2)

  if (is.numeric(Y)) {
    if (min(Y) < 0) Y <- Y - min(Y)
  } else {
    stopifnot(is.logical(Y))
  }

  stopifnot(
    setequal(levels(as.factor(A)), c(0, 1, 2)) |
      setequal(levels(as.factor(A)), c(0, 1))
  )

  owl_X <- X
  owl_Y <- (A > 0)
  owl_weights <- ifelse(is.matrix(propensity),
    as.numeric(Y) / propensity[cbind(1:length(Y), owl_Y + 1)],
    owl_weights <- as.numeric(Y) / propensity
  )
  aucs <- do.call(
    stabilityGLM,
    args = list(list(X = owl_X, Y = owl_Y, weights = owl_weights), owl_args)
  )

  return(aucs)
}

#' Implements the modified outcome approach
#'
#' @param A target variant
#' @param X rest of the genotype
#' @param Y phenotype
#' @param propensity propensity scores matrix
#' @param parallel whether to perform support estimation in a
#' parallelized fashion with the \code{\link{stabilityBIG}} function
#' @param ... additional arguments to be passed to \code{stabilityGLM} or
#' \code{stabilityBIG}
#'
#' @export
modified_outcome <- function(A, X, Y, propensity, parallel = FALSE, ...) {
  mod_args <- list(...)
  if ("family" %in% names(mod_args)) {
    stopifnot(mod_args["family"] == "gaussian")
  } else {
    mod_args["family"] <- "gaussian"
  }

  stopifnot(is.matrix(propensity))
  stopifnot(dim(propensity)[2] == 2)
  if (is.logical(Y)) Y <- 2 * as.numeric(Y) - 1

  mod_Y <- Y * ((A > 0) / propensity[, 2] - (A == 0) / propensity[, 1])

  if (parallel == TRUE) {
    mod_X <- X
    aucs <- do.call(
      stabilityBIG,
      args = list(list(X = mod_X, Y = mod_Y), mod_args)
    )
  } else {
    mod_X <- bigmemory::as.big.matrix(X, shared = FALSE)
    aucs <- do.call(
      stabilityGLM,
      args = list(list(X = mod_X, Y = mod_Y), mod_args)
    )
  }

  return(aucs)
}

#' Implements the normalized modified outcome approach
#'
#' @param A target variant
#' @param X rest of the genotype
#' @param Y phenotype
#' @param propensity propensity scores
#' @param parallel whether to perform support estimation in a
#' parallelized fashion
#' @param ... additional arguments to be passed to \code{stabilityGLM} or
#' \code{stabilityBIG}
#'
#' @export
normalized_outcome <- function(A, X, Y, propensity, parallel = FALSE, ...) {
  norm_args <- list(...)
  if ("family" %in% names(norm_args)) {
    stopifnot(norm_args["family"] == "gaussian")
  } else {
    norm_args["family"] <- "gaussian"
  }

  stopifnot(is.matrix(propensity))
  stopifnot(dim(propensity)[2] == 2)
  if (is.logical(Y)) Y <- 2 * as.numeric(Y) - 1

  propensity <- t(t(propensity) /
    c(
      sum(1 / (propensity[, 1][A == 0])),
      sum(1 / (propensity[, 2][A > 0]))
    ))

  norm_Y <- Y * ((A > 0) / propensity[, 2] -
    (A == 0) / propensity[, 1])

  if (parallel == TRUE) {
    norm_X <- X
    aucs <- do.call(
      stabilityBIG,
      args = list(list(X = norm_X, Y = norm_Y), norm_args)
    )
  } else {
    norm_X <- bigmemory::as.big.matrix(X, shared = FALSE)
    aucs <- do.call(
      stabilityGLM,
      args = list(list(X = norm_X, Y = norm_Y), norm_args)
    )
  }

  return(aucs)
}
#' Implements the shifted modified outcome approach
#'
#' @param A target variant
#' @param X rest of the genotype
#' @param Y phenotype
#' @param propensity propensity scores
#' @param parallel whether to perform support estimation in a
#' parallelized fashion
#' @param shift term to be added to the propensity scores to avoid
#' numerical stability
#' @param ... additional arguments to be passed to \code{stabilityGLM} or
#' \code{stabilityBIG}
#'
#' @export
shifted_outcome <- function(A, X, Y, propensity, shift = .1, parallel = FALSE, ...) {
  shift_args <- list(...)
  if ("family" %in% names(shift_args)) {
    stopifnot(shift_args["family"] == "gaussian")
  } else {
    shift_args["family"] <- "gaussian"
  }

  stopifnot(is.matrix(propensity))
  stopifnot(dim(propensity)[2] == 2)
  if (is.logical(Y)) Y <- 2 * as.numeric(Y) - 1

  shift_Y <- Y * ((A > 0) / (propensity[, 2] + shift) -
    (A == 0) / (propensity[, 1] + shift))

  if (parallel == TRUE) {
    shift_X <- X
    aucs <- do.call(
      stabilityBIG,
      args = list(list(X = shift_X, Y = shift_Y), shift_args)
    )
  } else {
    shift_X <- bigmemory::as.big.matrix(X, shared = FALSE)
    aucs <- do.call(
      stabilityGLM,
      args = list(list(X = shift_X, Y = shift_Y), shift_args)
    )
  }

  return(aucs)
}
#' Implements the robust modified outcome approach
#'
#' @param A target variant
#' @param X rest of the genotype
#' @param Y phenotype
#' @param propensity propensity scores
#' @param parallel whether to perform support estimation in a
#' parallelized fashion
#' @param ... additional arguments to be passed to \code{stabilityGLM} or
#' \code{stabilityBIG}
#'
#' @export
robust_outcome <- function(A, X, Y, propensity, parallel = FALSE, ...) {
  robust_args <- list(...)
  if ("family" %in% names(robust_args)) {
    stopifnot(robust_args["family"] == "gaussian")
  } else {
    robust_args["family"] <- "gaussian"
  }

  stopifnot(is.matrix(propensity))
  stopifnot(dim(propensity)[2] == 2)
  if (is.logical(Y)) Y <- 2 * as.numeric(Y) - 1

  C1 <- sum((as.numeric(A > 0) - propensity[, 2]) / propensity[, 2]) /
    sum(((as.numeric(A > 0) - propensity[, 2]) / propensity[, 2])^2)
  C0 <- sum((as.numeric(A == 0) - propensity[, 1]) / propensity[, 1]) /
    sum(((as.numeric(A == 0) - propensity[, 1]) / propensity[, 1])^2)

  propensity <- (propensity - matrix(rep(c(C0, C1), each = length(A)), ncol = 2)) /
    (propensity**2)

  propensity <- t(t(propensity) /
    c(
      sum(1 / (propensity[, 1][A == 0])),
      sum(1 / (propensity[, 2][A > 0]))
    ))

  robust_Y <- Y * ((A > 0) / propensity[, 2] -
    (A == 0) / propensity[, 1])

  if (parallel == TRUE) {
    robust_X <- X
    aucs <- do.call(
      stabilityBIG,
      args = list(list(X = robust_X, Y = robust_Y), robust_args)
    )
  } else {
    robust_X <- bigmemory::as.big.matrix(X, shared = FALSE)
    aucs <- do.call(
      stabilityGLM,
      args = list(list(X = robust_X, Y = robust_Y), robust_args)
    )
  }

  return(aucs)
}
