#' Samples causal SNPs with different effect types
#'
#' The sampled SNPs are combined in a list of character vectors with the
#' following names: target, marginal, inter1 and inter2. Through the
#' parameters \code{overlap_marg} and \code{overlap_inter}, the synergistic
#' SNPs with the target can have additional marginal/epistatic effects. The
#' SNPs are consecutively sampled in the following order: target, marginal,
#' inter1 and inter2. For each SNP, we iteratively sample until the picked
#' SNP candidate meets the constraints defined by \code{thresh_MAF} and
#' \code{window_size} (see Arguments for more details) or until the maximum
#' number of resamplings is reached. To avoid redundancy, we sample at most
#' one SNP per cluster.
#'
#' @section Warning:
#' Make sure to supply the SNP IDs in \code{names(clusters)}. The SNPs in
#' the output list are referenced by their names.
#'
#' @param nX number of SNPs interacting with the target variant
#' @param nY number of SNPs with marginal effects
#' @param nZ12 nummber of SNP pairs with epistatic effects
#' @param clusters vector of cluster memberships. Typically, the output
#'   of \code{\link[stats]{cutree}}. For ease of identification,
#'   the SNP IDs in \code{names(clusters)} are mandatory.
#' @param MAF vector of minor allele frequencies. The order of the SNPs in
#'   \code{MAF} is identical to that in \code{clusters}.
#' @param thresh_MAF lower-bound on the minor allele frequencies of
#'   causal SNPs. Rare variants are inherently difficult to recover. Assessing
#'   the retrieval performance on common variants better reflects the true
#'   performance of the epistasis detection algorithm.
#' @param window_size in number of clusters. Beside the target variant, the
#'   other causal SNPs are sampled outside of a window centered around the
#'   target. On each side of the target variant, the number of clusters to
#'   discard is \code{window_size}.
#' @param overlap_marg number of SNPs with both synergistic effects with the
#'   target and mariginal effects
#' @param overlap_inter number of SNPs with both synergistic effects with the
#'   target and additional epistatic effects
#' @param max_iter maximum number of resamplings for each SNP. If exceeded,
#'   the function generates an error
#'
#' @return list of character vectors corresponding to the causal SNP IDs. The
#'   output list entries are: target, marginal, inter1 and inter2. An
#'   epistatic pair is obtained from the combination of two SNPs with
#'   identical positions in \code{inter1} and \code{inter2}.
#'
#' @export
sample_SNP <- function(nX, nY, nZ12, clusters, MAF, thresh_MAF = 0.2,
                       window_size = 3, overlap_marg = 0, overlap_inter = 0, max_iter = 1e4) {
  stopifnot((thresh_MAF > 0) & (thresh_MAF <= 0.5))
  stopifnot(!is.null(names(clusters)))
  stopifnot(length(clusters) == length(MAF))
  stopifnot(overlap_marg + overlap_inter < nX)
  stopifnot(overlap_marg < nY)
  stopifnot(overlap_inter < nZ12)

  if (is.null(names(MAF))) {
    names(MAF) <- names(clusters)
  } else {
    stopifnot(all(names(MAF) == names(clusters)))
  }
  active <- rep(0, 1 + nX + nY + 2 * nZ12 - (overlap_marg + overlap_inter))
  for (i in 1:length(active)) {
    iter <- 0
    cdt <- FALSE
    while (cdt == FALSE) {
      iter <- iter + 1
      candidate_SNP <- sample.int(length(clusters), size = 1)
      candidate_name <- names(clusters)[candidate_SNP]
      cdt <- (MAF[candidate_SNP] > thresh_MAF)
      if (i == 1) {
        cluster_target <- clusters[candidate_name]
      } else {
        cdt <- cdt & (abs(clusters[candidate_name] - cluster_target) > window_size)
      }
      if (iter > max_iter) stop("The MAF constraint can not be satisfied")
    }
    active[i] <- candidate_name
    MAF <- subset(MAF, (clusters != clusters[candidate_SNP]))
    clusters <- subset(clusters, (clusters != clusters[candidate_SNP]))
  }

  marg_idx <- sample.int(nX, overlap_marg)
  if (overlap_marg == 0) {
    inter_idx <- sample((1:nX), overlap_inter)
  }
  else {
    inter_idx <- sample((1:nX)[-marg_idx], overlap_inter)
  }
  SNP_list <- list(
    target = active[1], syner = active[2:(nX + 1)],
    marginal = active[c((nX + 2):(nX + nY + 1 - overlap_marg), 1 + marg_idx)],
    inter1 = active[c((nX + nY + 2):(nX + nY + nZ12 + 1 - overlap_inter) -
      overlap_marg, 1 + inter_idx)],
    inter2 = active[((nX + nY + nZ12 + 2 - (overlap_marg + overlap_inter)):
    length(active))]
  )

  return(SNP_list)
}

#' Samples effect sizes for the disease model
#'
#' The generated disease model is the list of effect size coeffcients.
#' The list comprises the following fields: "syner", "marg" and "inter".
#' "syner" is itself a list of numeric vectors with two entries named
#' "A0" and "A1". "A0" refers to the vector of effect sizes when the target
#' variant \eqn{A = 0}{A = 0}. Similarly, "A1" refers to the vector of effect
#' sizes in the case \eqn{A = 1}{A = 1}. The two other entries "marg" and
#' "inter" are, respectively, the marginal and epistatic effect sizes. The
#' effect sizes are independent and normally-distributed. The \code{mean}
#' parameter is either a list of vectors or a vector of length 4. If
#' \code{mean} is a vector, then the effet sizes for each type of effects are
#' identically distributed. Otherwise, the corresponding vector in the list
#' specifies their individual means. The same logic applies to \code{sd}, the
#' standard deviation parameter. For coherence, the parameters \code{mean} and
#' \code{sd} are encoded in the same order as the output.
#'
#' @param nX number of SNPs interacting with the target variant
#' @param nY number of SNPs with marginal effects
#' @param nZ12 nummber of SNP pairs with epistatic effects
#' @param mean vector or list of means
#' @param sd vector or list of standard deviations
#'
#' @return a list of vectors corresponding to the effect size coefficients.
#'
#' @export
gen_model <- function(nX, nY, nZ12, mean = rep(0, 4), sd = rep(1, 4)) {
  stopifnot(length(mean) == 4)
  stopifnot(length(sd) == 4)
  if (is.list(mean)) {
    stopifnot(
      (length(mean[[1]]) == nX) | (length(mean[[2]]) == nX) |
        (length(mean[[3]]) == nY) | (length(mean[[4]]) == nZ12)
    )
  }
  if (is.list(sd)) {
    stopifnot(
      (length(sd[[1]]) == nX) | (length(sd[[2]]) == nX) |
        (length(sd[[3]]) == nY) | (length(sd[[4]]) == nZ12)
    )
  }

  model <- list(
    syner = list(
      A0 = stats::rnorm(nX, mean = mean[[1]], sd = sd[[1]]),
      A1 = stats::rnorm(nX, mean = mean[[2]], sd = sd[[2]])
    ),
    marg = stats::rnorm(nY, mean = mean[[3]], sd = sd[[3]]),
    inter = stats::rnorm(nZ12, mean = mean[[4]], sd = sd[[4]])
  )

  return(model)
}

#' Simulates a binary phenotye
#'
#' The phenotypes are simulated according to a logistic regression model.
#' Depending on the chosen configuration in \code{\link{sample_SNP}}, the model
#' includes different effect types: synergistic effects with the target,
#' marginal effects and additional epistatic effects. As the generated
#' phenotypes vector is often to supplied to supervized learning tasks, we
#' offer the option to obtain a balanced dataset, through the \code{intercept}
#' paramter
#'
#' @param X genotype matrix
#' @param causal causal SNPs.
#' @param model disease model
#' @param intercept binary flag. If \code{intercept=TRUE}, a non-null intercept
#'   is added so that the output is (approximately) balanced between cases and
#'   controls.
#'
#' @return A vector of simulated phenotypes which are encoded as a two-level
#' factor (TRUE/FALSE).
#'
#' @export
sim_phenotype <- function(X, causal, model, intercept = TRUE) {
  stopifnot(!is.null(colnames(X)))
  stopifnot(!(is.null(causal[["target"]]) | is.null(causal[["syner"]]) |
    is.null(causal[["marginal"]]) | is.null(causal[["inter1"]]) |
    is.null(causal[["inter2"]])
  ))
  stopifnot(!(is.null(model[["syner"]][["A0"]]) |
    is.null(model[["syner"]][["A1"]]) |
    is.null(model[["marg"]]) | is.null(model[["inter"]])
  ))
  stopifnot(
    (length(causal[["syner"]]) == length(model[["syner"]][["A0"]])) &
      (length(causal[["syner"]]) == length(model[["syner"]][["A1"]])) &
      (length(causal[["marginal"]]) == length(model[["marg"]])) &
      (length(causal[["inter1"]]) == length(model[["inter"]])) &
      (length(causal[["inter2"]]) == length(model[["inter"]]))
  )

  risk <- (X[, causal$target] == 0) * (X[, causal$syner] %*% model$syner$A0) +
    (X[, causal$target] > 0) * (X[, causal$syner] %*% model$syner$A1) +
    (X[, causal$marginal] %*% model$marg) +
    (X[, causal$inter1] * X[, causal$inter2]) %*% model$inter

  if (intercept) {
    risk <- risk - mean(risk)
  }

  phenotypes <- stats::runif(dim(X)[1]) < (1 / (1 + exp(-risk)))

  return(phenotypes)
}

#' Merges a number of clusters around the target
#'
#' Replaces the indices of neighbor clusters with \code{center}, the
#' target cluster index. The neighborhood is defined according to the parameter
#' \code{k} (see Arguments for more details). The real purpose of the function
#' \code{\link{merge_cluster}} is to define an enlarged window of SNPs which
#' are in linkage disequilibrium with the target. Subsequently, we filter
#' them out for the estimate of the propensity scores.
#'
#' @param clusters vector of cluster memberships. Typically, the output
#'   of \code{\link[stats]{cutree}}
#' @param center the target variant cluster
#' @param k vector or integer. if \code{k} is a vector, it corresponds to
#'   the cluster indices to be updated. Otherwise, \code{k} is an integer
#'   and the cluster indices to be updated lie between \code{center-k} and
#'   \code{center+k}.
#'
#' @return The updated cluster membership vector. The cluster indexing is also
#'   updated so that the maximum cluster index is equal to the total number of
#'   clusters after merging.
#'
#' @examples
#' hc <- hclust(dist(USArrests))
#' clusters <- cutree(hc, k = 10)
#' merge_cluster(clusters, center=5, k=2)
#'
#' @export
merge_cluster <- function(clusters, center, k = 3) {
  stopifnot(center %in% clusters)
  if (is.atomic(k) & length(k) == 1L) {
    stopifnot((2 * k + 1) <= nlevels(as.factor(clusters)))
    idx_window <- center - (-k:k)
  } else {
    stopifnot(all(k %in% clusters))
    idx_window <- k
  }

  clusters[which(clusters %in% idx_window)] <- center
  clusters <- vapply(clusters, function(s) match(s, sort(unique(clusters))))

  return(clusters)
}
