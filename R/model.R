#' Samples causal SNPs
#'
#'
#'
#'
#' @param nX
#' @param nY
#' @param nZ12
#' @param clusters
#' @param MAF
#' @param thresh.MAF
#' @windo
#'
#'
#'
#'
#'
#'
sample_SNP <- function(nX, nY, nZ12, clusters, MAF, thresh.MAF = 0.2,
                       window.size = 3, overlap.marg = 0, overlap.inter = 0) {
  stopifnot(overlap.marg + overlap.inter < nX)
  active <- rep(NULL, 1 + nX + nY + 2 * nZ12 - (overlap.marg + overlap.inter))

  copy.clusters <- clusters
  copy.MAF <- MAF

  max.iter <- 1e4

  for (i in 1:(1 + nX + nY + 2 * nZ12 - overlap.inter - overlap.marg)) {
    iter <- 0
    cdt <- FALSE
    while (cdt == FALSE) {
      iter <- iter + 1
      candidate.SNP <- sample.int(length(copy.clusters), size = 1)
      candidate.name <- names(copy.clusters)[candidate.SNP]
      cdt <- (copy.MAF[candidate.SNP] > thresh.MAF)
      if (i > 1) cdt <- cdt & (abs(clusters[candidate.name] - clusters[active[1]]) > window.size)
      if (iter > max.iter) stop("The MAF constraint can not be satisfied")
    }
    active[i] <- candidate.name
    copy.MAF <- copy.MAF[copy.clusters != copy.clusters[candidate.SNP]]
    copy.clusters <- copy.clusters[copy.clusters != copy.clusters[candidate.SNP]]
  }


  marg.idx <- sample.int(nX, overlap.marg)
  if (overlap.marg == 0) {
    inter.idx <- sample((1:nX), overlap.inter)
  }
  else {
    inter.idx <- sample((1:nX)[-marg.idx], overlap.inter)
  }

  return(list(
    target = active[1], syner = active[2:(nX + 1)],
    marginal = active[c((nX + 2):(nX + nY + 1 - overlap.marg), 1 + marg.idx)],
    inter1 = active[c((nX + nY + 2):(nX + nY + nZ12 + 1 - overlap.inter) -
                        overlap.marg, 1 + inter.idx)],
    inter2 = active[((nX + nY + nZ12 + 2 - (overlap.marg + overlap.inter)):
                       length(active))]
  ))
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
  stopifnot(length(mean)==4)
  stopifnot(length(sd)==4)
  if(is.list(mean)){
    stopifnot(
      (length(mean[[1]])==nX)|(length(mean[[2]])==nX)|
        (length(mean[[3]])==nY)|(length(mean[[4]])==nZ12)
    )
  }
  if(is.list(sd)){
    stopifnot(
      (length(sd[[1]])==nX)|(length(sd[[2]])==nX)|
        (length(sd[[3]])==nY)|(length(sd[[4]])==nZ12)
    )
  }

  model <- list(
    syner = list(
      A0 = rnorm(nX, mean = mean[[1]], sd = sd[[1]]),
      A1 = rnorm(nX, mean = mean[[2]], sd = sd[[2]])
    ),
    marg = rnorm(nY, mean = mean[[3]], sd = sd[[3]]),
    inter = rnorm(nZ12, mean = mean[[4]], sd = sd[[4]])
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
sim_phenotype <- function(X, causal, model, intercept=TRUE) {
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

  if(intercept){
    risk <- risk - mean(risk)
  }

  phenotypes <- runif(dim(X)[1]) < (1 / (1 + exp(-risk)))

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
#'   of \code{\link[stats::cutree]{cutree}}
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
#' merge.cluster(clusters, center=5, window.size=2)
#'
#' @export
merge_cluster <- function(clusters, center, k = 3) {
  stopifnot(!is.integer(clusters))
  stopifnot(center %in% clusters)
  if (is.atomic(k)) {
    stopifnot((2 * k + 1) <= nlevels(as.factor(clusters)))
    idx_window <- center - (-k:k)
  } else {
    stopifnot(all(k %in% clusters))
    idx_window <- k
  }

  clusters[which(clusters %in% idx_window)] <- center
  clusters <- sapply(clusters, function(s) match(s, sort(unique(clusters))))

  return(clusters)
}
