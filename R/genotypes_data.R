#' Arabidopsis QTL data on gravitropism
#'
#' Data from a QTL experiment on gravitropism in
#' Arabidopsis, with data on 162 recombinant inbred lines (Ler x
#' Cvi). The outcome is the root tip angle (in degrees) at two-minute
#' increments over eight hours.
#'
#' @docType data
#'
#' @usage data(genotypes)
#'
#' @name genotypes
#'
#' @format An object of class \code{"cross"}; see.
#'
#' @keywords datasets
#'
#' @references Su, Z., Marchini, J., & Donnelly, P. (2011). HAPGEN2: Simulation
#' of multiple disease SNPs. Bioinformatics, 27(16), 2304–2305.
#'
#' @examples
#' data(genotypes)
#' \dontrun{fast_HMM(genotypes, fp_path = "~/path/to/fastPHASE")}
"genotypes"
