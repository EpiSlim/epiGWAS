#' Simulated genotypes
#'
#' We simulated 300 unphased european genotypes. For that matter, we
#' used the HAPGEN2 software and the 1000 genome phase 3 reference data.
#' The simulated regions are located on the 22 chromosome between the
#' nucleotide positions 16061016 (rs9617528) and 19976834 (rs887201).
#' Only the markers of the Affymetrix 500K are included.
#'
#' @docType data
#'
#' @usage data(genotypes)
#'
#' @name genotypes
#'
#' @format An integer matrix with 300 rows and 450 columns. The SNP rsIDs
#' and positions, in addition to their reference and alternate alleles, are
#' combined in \code{colnames(X)}.
#'
#' @keywords datasets
#'
#' @references Su, Z., Marchini, J., & Donnelly, P. (2011). HAPGEN2: Simulation
#' of multiple disease SNPs. Bioinformatics, 27(16), 2304–2305.
#'
#' @references Consortium, T. 1000 G. P., Auton, A., Abecasis, G. R., Altshuler
#' (Co-Chair), D. M., Durbin (Co-Chair), R. M., Abecasis, G. R., … Abecasis,
#' G. R. (2015). A global reference for human genetic variation. Nature, 526, 68.
#'
#' @examples
#' data(genotypes)
#' \dontrun{hmm <- fast_HMM(genotypes, fp_path = '~/path/to/fastPHASE')}
"genotypes"
