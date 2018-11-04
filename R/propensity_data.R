#' propensity scores
#'
#' All SNPs in the genomic interval [18038910, 18288361] are removed
#' from \code{genotypes} to alleviate linkage disequilibrium around the
#' target rs2535708 (position 18184169). Afterwards, we apply \code{\link{fast_HMM}}
#' and then \code{\link{cond_prob}} to obtain the propensity scores. The results
#' are stored in \code{propensity}. 
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
4
4
