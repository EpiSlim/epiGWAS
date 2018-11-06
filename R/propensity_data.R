#' propensity scores
#'
#' To obtain the scores for the samples in \code{\link{genotypes}},
#' we first remove all SNPs in the genomic interval
#' [18038910, 18288361] to alleviate linkage disequilibrium around the
#' target SNP rs2535708 (position 18184169). Afterwards, we apply
#' \code{\link{fast_HMM}} followed by \code{\link{cond_prob}}.
#' The results are stored in \code{propensity}.
#'
#' @docType data
#'
#' @usage data(propensity)
#'
#' @name propensity
#'
#' @format a numeric vector
#'
#' @examples
#' data(propensity)
"propensity"
