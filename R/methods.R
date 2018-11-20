#' Implements the outcome weighted learning approach
#'
#' To recover the synergistic interactions between the target \code{A} and the
#' rest of the genotype \code{X}, \code{OWL} formulates a weighted binary
#' classification problem. The outcome is the mapping of \code{A} to \{0,1\}. The
#' covariates are \code{X}. The propensity scores and the phenotypes are
#' combined in the sample weights \eqn{Y/\pi(A\lvert X)}{Y/P(A|X)}. For binary
#' phenotypes, OWL is a case-only approach. The approach also accommodates
#' nonnegative continuous phenotypes.
#'
#' @param A target variant. If not binary, the variable A must be encoded
#' as either (0, 1) or (0, 1, 2).
#' @param X rest of the genotype
#' @param Y phenotype (binary or continuous)
#' @param propensity propensity scores (a vector or a two-column matrix)
#' @param ... additional arguments to \code{\link{stabilityGLM}}
#'
#' @return a vector containing the area under the stability selection path for
#'   each variable in \code{X}
#'
#' @details For continuous phenotypes, if the outcome \code{Y} is
#'   not nonnegative, it is translated to make it nonnegative.
#'
#' @references Zhao, Y., Zeng, D., Rush, A. J., & Kosorok, M. R. (2012).
#'   Estimating Individualized Treatment Rules Using Outcome Weighted Learning.
#'   Journal of the American Statistical Association, 107(499), 1106–1118.
#'
#' @export
OWL <- function(A, X, Y, propensity, ...) {
  owl_args <- list(...)
  if ("family" %in% names(owl_args)) {
    stopifnot(owl_args["family"] == "binomial")
  } else {
    owl_args["family"] <- "binomial"
  }

  stopifnot(all(propensity >= 0) & (propensity <= 1))
  if (is.matrix(propensity)) {
    stopifnot(dim(propensity)[2] == 2)
  }

  if (is.numeric(Y)) {
    if (min(Y) < 0) {
      Y <- Y - min(Y)
    }
  } else {
    stopifnot(is.logical(Y))
  }

  stopifnot(setequal(levels(as.factor(A)), c(0, 1, 2)) | setequal(
    levels(as.factor(A)),
    c(0, 1)
  ) | is.logical(A))

  owl_X <- X
  if (is.logical(A)) {
    owl_Y <- A
  } else {
    owl_Y <- (A > 0)
  }
  if (is.matrix(propensity)) {
    owl_weights <- as.numeric(Y) / propensity[cbind(
      seq_len(length(Y)),
      owl_Y + 1
    )]
  } else {
    owl_weights <- as.numeric(Y) / propensity
  }
  aucs <- do.call(stabilityGLM, args = append(list(
    X = owl_X, Y = owl_Y,
    weights = owl_weights
  ), owl_args))

  return(aucs)
}

#' Implements the modified outcome approach
#'
#' In the modified outcome approach, we estimate the risk difference
#' \eqn{\mathbb{E}\left[Y\lvert A=1,X\right]-\mathbb{E}\left[Y\lvert A=0,X\right]}{E[Y|A=1,X]-E[Y|A=0,X]}.
#' The risk difference measures the synergy between \code{A} and the set of
#' covariates in \code{X}. For genome-wide association studies, it can be
#' interpreted as a pure epistatic term. However, for a single sample, we only
#' observe one of the two possibilities A=1 or A=0, making the direct
#' estimate of the risk difference impossible. Through propensity scores,
#' modified outcome was proposed as a solution to this problem. The risk
#' difference is recovered by constructing a modified outcome that combines
#' A, Y and the propensity score \eqn{\pi(A\lvert X)}{P(A|X)}:
#' \eqn{Y \times \left[\frac{A}{\pi(A=1\lvert X)} - \frac{1-A}{1-\pi(A=1\lvert X)} \right]}{Y x [A/P(A=1|X) -(1-A)/P(A=0|X)]}.
#' The use of \code{\link{stabilityGLM}} or \code{\link{stabilityBIG}} for
#' the modified outcome regression allows us to recover the
#' interacting components within \code{X}.
#'
#' @param A target variant
#' @param X rest of the genotype
#' @param Y phenotype
#' @param propensity propensity scores vector/matrix. If given as a matrix,
#' the first column is \eqn{\pi(A = 0\lvert X)}{P(A = 0|X)} while the second
#' is \eqn{\pi(A = 1\lvert X)}{P(A = 1|X)}
#' @param parallel whether to perform support estimation in a
#' parallelized fashion with the \code{\link{stabilityBIG}} function
#' @param ... additional arguments to be passed to \code{stabilityGLM} or
#' \code{stabilityBIG}
#'
#' @return a vector containing the area under the stability selection path for
#'   each variable in \code{X}
#'
#' @references Rosenbaum, Paul R., and Donald B. Rubin. 'The central role of
#' the propensity score in observational studies for causal effects.'
#' Biometrika 70.1 (1983): 41-55.
#'
#' @export
modified_outcome <- function(A, X, Y, propensity, parallel = FALSE, ...) {
  mod_args <- list(...)
  if ("family" %in% names(mod_args)) {
    stopifnot(mod_args["family"] == "gaussian")
  } else {
    mod_args["family"] <- "gaussian"
  }

  if (is.vector(propensity)) {
    propensity <- cbind(
      (A == 0) * propensity + (A > 0) * (1 - propensity),
      (A > 0) * propensity + (A == 0) * (1 - propensity)
    )
  }
  if (is.logical(Y)) {
    Y <- 2 * as.numeric(Y) - 1
  }

  mod_Y <- Y * ((A > 0) / propensity[, 2] - (A == 0) / propensity[, 1])

  if (parallel == TRUE) {
    mod_X <- bigmemory::as.big.matrix(X, shared = FALSE, type = "double")
    aucs <- do.call(stabilityBIG, args = append(
      list(X = mod_X, Y = mod_Y),
      mod_args
    ))
  } else {
    mod_X <- X
    aucs <- do.call(stabilityGLM, args = append(
      list(X = mod_X, Y = mod_Y),
      mod_args
    ))
  }

  return(aucs)
}

#' Implements the normalized modified outcome approach
#'
#' Normalized modified outcome is an improvement to \code{\link{modified_outcome}}.
#' Its large-sample variance is lower than the original modified outcome approach.
#' The only difference between the two methods lies in the normalization of the
#' propensity scores. The inverses of the propensity scores
#' \eqn{1/\pi(A=1\lvert X)}{1/P(A=1|X)} and \eqn{1/\pi(A=0\lvert X)}{1/P(A=0|X)} are
#' respectively normalized by their sum
#' \eqn{\sum_{i} 1/\pi(A_i=1\lvert X_i)}{sum _i 1/P(A_i=1|X_i)} and
#' \eqn{\sum_{i} 1/\pi(A_i=0\lvert X_i)}{sum _i 1/P(A_i=0|X_i)}.
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
#' @return a vector containing the area under the stability selection path for
#'   each variable in \code{X}
#'
#' @export
normalized_outcome <- function(A, X, Y, propensity, parallel = FALSE, ...) {
  norm_args <- list(...)
  if ("family" %in% names(norm_args)) {
    stopifnot(norm_args["family"] == "gaussian")
  } else {
    norm_args["family"] <- "gaussian"
  }

  if (is.vector(propensity)) {
    propensity <- cbind(
      (A == 0) * propensity + (A > 0) * (1 - propensity),
      (A > 0) * propensity + (A == 0) * (1 - propensity)
    )
  }
  if (is.logical(Y)) {
    Y <- 2 * as.numeric(Y) - 1
  }

  propensity <- t(t(propensity) / c(sum(1 / (propensity[, 1][A == 0])), sum(1 / (propensity[
    ,
    2
  ][A > 0]))))

  norm_Y <- Y * ((A > 0) / propensity[, 2] - (A == 0) / propensity[, 1])

  if (parallel == TRUE) {
    norm_X <- bigmemory::as.big.matrix(X, shared = FALSE, type = "double")
    aucs <- do.call(stabilityBIG, args = append(
      list(X = norm_X, Y = norm_Y),
      norm_args
    ))
  } else {
    norm_X <- X
    aucs <- do.call(stabilityGLM, args = append(
      list(X = norm_X, Y = norm_Y),
      norm_args
    ))
  }

  return(aucs)
}
#' Implements the shifted modified outcome approach
#'
#' Shifted modified outcome is an improvement to modified outcome. It is
#' a heuristic which consists in the addition of of a small translation term to
#' the inverse of the propensity score. The goal is to avoid numerical
#' instability due to low propensity scores values. More precisely, the
#' inverses of the propensity scores become
#' \eqn{1/(\pi(A\lvert X) + \xi)}{1/(P(A|X) + shift)}. We recommend keeping
#' the default value of the parameter \eqn{\xi}{shift} at 0.1.
#'
#' @param A target variant
#' @param X rest of the genotype
#' @param Y phenotype
#' @param propensity propensity scores
#' @param parallel whether to perform support estimation in a
#' parallelized fashion
#' @param shift regularization term to be added to the propensity scores to
#' avoid numerical stability
#' @param ... additional arguments to be passed to \code{stabilityGLM} or
#' \code{stabilityBIG}
#'
#' @return a vector containing the area under the stability selection path for
#'   each variable in \code{X}
#'
#' @export
shifted_outcome <- function(A, X, Y, propensity, parallel = FALSE, shift = 0.1,
                            ...) {
  shift_args <- list(...)
  if ("family" %in% names(shift_args)) {
    stopifnot(shift_args["family"] == "gaussian")
  } else {
    shift_args["family"] <- "gaussian"
  }

  if (is.vector(propensity)) {
    propensity <- cbind(
      (A == 0) * propensity + (A > 0) * (1 - propensity),
      (A > 0) * propensity + (A == 0) * (1 - propensity)
    )
  }
  if (is.logical(Y)) {
    Y <- 2 * as.numeric(Y) - 1
  }

  shift_Y <- Y * ((A > 0) / (propensity[, 2] + shift) - (A == 0) / (propensity[
    ,
    1
  ] + shift))

  if (parallel == TRUE) {
    shift_X <- bigmemory::as.big.matrix(X, shared = FALSE, type = "double")
    aucs <- do.call(stabilityBIG, args = append(
      list(X = shift_X, Y = shift_Y),
      shift_args
    ))
  } else {
    shift_X <- X
    aucs <- do.call(stabilityGLM, args = append(
      list(X = shift_X, Y = shift_Y),
      shift_args
    ))
  }

  return(aucs)
}
#' Implements the robust modified outcome approach
#'
#' A key feature of \code{robust_outcome} is its resilience to the misspecification
#' of propensity scores, which is a major limitation of classical
#' modified outcome approaches.
#' Except for \code{shifted_outcome}, all of
#' the modified outcome approaches belong to a parameterized class of unbiased
#' estimators for the risk difference term
#' \eqn{\mathbb{E}\left[Y\lvert A=1,X\right]-\mathbb{E}\left[Y\lvert A=0,X\right]}{E[Y|A=1,X]-E[Y|A=0,X]}.
#' Within that class, robust modified outcome is the approach with the least
#' large-sample variance. For more details about this approach, see Lunceford
#' and Davidian (2004)
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
#' @return a vector containing the area under the stability selection path for
#'   each variable in \code{X}
#'
#' @references Lunceford, J. K., & Davidian, M. (2004). Stratification and
#' weighting via the propensity score in estimation of causal treatment effects:
#' A comparative study. Statistics in Medicine, 23(19), 2937–2960.
#'
#' @export
robust_outcome <- function(A, X, Y, propensity, parallel = FALSE, ...) {
  robust_args <- list(...)
  if ("family" %in% names(robust_args)) {
    stopifnot(robust_args["family"] == "gaussian")
  } else {
    robust_args["family"] <- "gaussian"
  }

  if (is.vector(propensity)) {
    propensity <- cbind(
      (A == 0) * propensity + (A > 0) * (1 - propensity),
      (A > 0) * propensity + (A == 0) * (1 - propensity)
    )
  }
  if (is.logical(Y)) {
    Y <- 2 * as.numeric(Y) - 1
  }

  C1 <- sum((as.numeric(A > 0) - propensity[, 2]) / propensity[, 2]) / sum(((as.numeric(A >
    0) - propensity[, 2]) / propensity[, 2])^2)
  C0 <- sum((as.numeric(A == 0) - propensity[, 1]) / propensity[, 1]) / sum(((as.numeric(A ==
    0) - propensity[, 1]) / propensity[, 1])^2)

  propensity <- (propensity - matrix(rep(c(C0, C1), each = length(A)),
    ncol = 2
  )) / (propensity^2)

  propensity <- t(t(propensity) / c(sum(1 / (propensity[, 1][A == 0])), sum(1 / (propensity[
    ,
    2
  ][A > 0]))))

  robust_Y <- Y * ((A > 0) / propensity[, 2] - (A == 0) / propensity[, 1])

  if (parallel == TRUE) {
    robust_X <- bigmemory::as.big.matrix(X, shared = FALSE, type = "double")
    aucs <- do.call(stabilityBIG, args = append(
      list(X = robust_X, Y = robust_Y),
      robust_args
    ))
  } else {
    robust_X <- X
    aucs <- do.call(stabilityGLM, args = append(
      list(X = robust_X, Y = robust_Y),
      robust_args
    ))
  }

  return(aucs)
}


#' Runs a selection of epistasis detection methods in a joint manner
#'
#' This function is a wrapper for the different epistasis detection
#' methods implemented in this package. If \code{methods} is \code{"all"}, we run
#' \code{OWL} and the four modified outcome approaches. Otherwise, we run a
#' selection of those methods. In this case, the \code{methods} argument is a
#' character vector with its entries being the names of the functions to call.
#'
#' @param A target variant
#' @param X rest of the genotype
#' @param Y phenotype
#' @param propensity propensity scores
#' @param parallel whether to perform support estimation in a
#'     parallelized fashion for the modified outcome family of methods
#' @param methods character vector for the epistasis detection methods to call
#' @param shift regularization parameter for \code{\link{shifted_outcome}}
#' @param ... additional arguments to be passed to \code{stabilityGLM} or
#'     \code{stabilityBIG}
#'
#' @return list of numeric vectors. Each vector corresponds to the auc scores
#'     of a particular method in \code{methods}.
#'
#' @seealso \code{OWL}, \code{modified_outcome}, \code{shifted_outcome},
#'     \code{normalized_outcome} and \code{robust_outcome}
#'
#' @export
epiGWAS <- function(A, X, Y, propensity, methods = "all", parallel = TRUE, shift = 0.1, ...) {
  methods_full <- c(
    "OWL", "modified_outcome", "normalized_outcome",
    "shifted_outcome", "robust_outcome"
  )
  if (any(methods == "all")) {
    methods <- methods_full
  } else {
    stopifnot(all(methods %in% methods_full))
  }

  aucs <- list()
  for (m in methods) {
    method_args <- append(list(A = A, X = X, Y = Y, propensity = propensity), list(...))
    if (m == "shifted_outcome") method_args <- append(method_args, list(shift = shift))
    if (m != "OWL") {
      method_args <- append(method_args, list(parallel = parallel))
    } else {
      if ("ncores" %in% names(method_args)) method_args[["ncores"]] <- NULL
    }
    aucs[[m]] <- do.call(m, args = method_args)
  }
  return(aucs)
}
