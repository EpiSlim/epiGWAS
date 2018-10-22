#' Implements BOOST SNP-SNP interaction test
#'
#' For a pair of SNPs (\eqn{X_1}{X1}, \eqn{X_2}{X2}) and a binary phenotype
#' \eqn{Y}{Y}, the \code{\link{BOOST}} function computes the ratio of maximum
#' log-likelihoods for two models: the full model and the main effect model.
#' Mathematically speaking, the full model is a logistic regression model with
#' both main effect terms and interaction terms
#' \eqn{\left(X_1, X_2, X_1\times X_2\right)} {(X1, X2, X1, X1 x X2)}.
#' The main effect model is a logistic regression model with only
#' \eqn{\left(X_1, X_2\right)} {(X1, X2)} as covariates. Since we are
#' interested in the synergies with a single variant, we do not implement
#' the initial sure screening stage which filters out non-significant pairs.
#'
#' @seealso the webpage \url{http://bioinformatics.ust.hk/BOOST.html} provides.
#' additional details about the BOOST software
#'
#' @param A target variant. The SNP A is encoded as 0, 1, 2.
#' @param X genotype matrix (excluding A). The only accepted SNP values are
#' also 0, 1 and 2.
#' @param Y observed phenotype. A two-level factor.
#' @param ncores number of threads (default 1)
#'
#' @return the interaction statistic for each column in \code{X} and \code{A}
#'
#' @examples
#' X <- matrix((runif(5e2, min = 0, max = 1) < 0.5) + (runif(5e2, min = 0, max = 1) < 0.5), nrow = 50)
#' A <- (runif(50, min = 0, max = 1) < 0.5) + (runif(50, min = 0, max = 1) < 0.5)
#' Y <- runif(50, min = 0, max = 1) < 1/(1+exp(-.5 * A * X[, 3] + .25 * A * X[, 7]))
#' BOOST(A, X, Y)
#'
#' @export
BOOST <- function(A, X, Y, ncores = 1) {
  stopifnot(ncores <= parallel::detectCores())
  stopifnot(dim(X)[1] == length(A))
  stopifnot(dim(X)[1] == length(Y))
  stopifnot(setequal(levels(as.factor(A)), c(0, 1, 2)))
  stopifnot(setequal(levels(as.factor(X)), c(0, 1, 2)))
  stopifnot(is.factor(Y))
  stopifnot(length(levels(Y)) == 2)

  # Internal function for computing the likelihood ratio statistic
  ratio <- function(z) {
    data_AX <- data.frame(x1 = factor(A), x2 = factor(z), y = Y)
    fit01 <- glm(y ~ x1 + x2, family = "binomial", data = data_AX)
    fit2 <- glm(y ~ x1 + x2 + x1 * x2, family = "binomial", data = data_AX)

    return(fit01$dev - fit2$dev)
  }

  if (ncores == 1) {
    loglikelihood <- apply(X, 2, ratio)

  } else {
    cl <- snow::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)

    loglikelihood <-parallel::parApply(cl, X, 2, ratio)

    snow::stopCluster()

  }

  return(loglikelihood)
}
