#' SNP minor allele frequencies
#'
#' \code{maf} contains the minor allele frequencies for the 450 SNPs in
#' the \code{\link{genotypes}} dataset
#'
#' @docType data
#'
#' @usage data(maf)
#'
#' @name maf
#'
#' @format a numeric vector
#'
#' @examples
#' data(maf)
#' all((maf <= 0.5) & (maf >= 0))
"maf"
