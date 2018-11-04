#' fastPHASE hidden Markov model for the genotypes dataset
#'
#' Since the fastPHASE executable is not part of this package, we instead
#' provide its output in order to demonstrate the full features of the
#' epiGWAS package in the accompanying vignette. By means of the wrapper
#' function \code{\link{fast_HMM}},  we ran fastPHASE with the
#' paramaterization \code{n_state = 10} and  \code{n_iter = 25}. The
#' structure of the \code{hmm} list is explained in greater detail in
#' the documentation of \code{\link[SNPknock]{SNPknock.fp.loadFit_hmm}}
#'
#' @docType data
#'
#' @usage data(hmm)
#'
#' @name hmm
#'
#' @format A list consisting of three fields: \code{pInit}, \code{Q} and
#' \code{pEmit}
#'
#' @keywords datasets
#'
#' @references Scheet, P., & Stephens, M. (2006). A fast and flexible
#'     statistical model for large-scale population genotype data: applications
#'     to inferring missing genotypes and haplotypic phase. American Journal of
#'     Human Genetics, 78(4), 629–644. https://doi.org/10.1086/502802
#'
#' @seealso \code{\link{fast_HMM}} and \code{\link[SNPknock]{SNPknock.fp.loadFit_hmm}}
#'
#' @examples
#' data(hmm)
#' setequal(names(hmm), c("pEmit", "Q", "pInit"))
"hmm"
4
4
