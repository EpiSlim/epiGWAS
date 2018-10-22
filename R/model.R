sample.SNP <- function(nX, nY, nZ12, clusters, MAF, thresh.MAF=0.2, window.size =3, overlap.marg=0, overlap.inter=0) {
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
    inter1 = active[c((nX + nY + 2):(nX + nY + nZ12 + 1 - overlap.inter) - overlap.marg, 1 + inter.idx)],
    inter2 = active[((nX + nY + nZ12 + 2 - (overlap.marg + overlap.inter)):length(active))]
  ))
}

gen.model <- function(nX, nY, nZ12, mode=1, mean=rep(0, 4), sd=rep(1, 4)) {
  stopifnot((length(mean) - mode == 3) & (length(sd) - mode == 3))
  return(list(
    syner = list(A0 = rnorm(nX, mean = mean[1], sd = sd[1]), A1 = rnorm(nX, mean = mean[2], sd = sd[2])),
    marg = rnorm(nY, mean = mean[3], sd = sd[3]),
    inter = rnorm(nZ12, mean = mean[4], sd = sd[4]),
    direct = switch(mode, 0, rnorm(1, mean = mean[5], sd = sd[5])),
    mode = c("Synergistic only", "Combined effects")[mode]
  ))
}

sim.phen <- function(X, active, model) {
  risk <- (X[, active$target] == 0) * (X[, active$syner] %*% model$syner$A0) +
    (X[, active$target] > 0) * (X[, active$syner] %*% model$syner$A1) +
    (X[, active$marginal] %*% model$marg) + (X[, active$inter1] * X[, active$inter2]) %*% model$inter +
    (X[, active$target] > 0) * model$direct

  risk <- risk - mean(risk)
  return(runif(nrow(X)) < (1 / (1 + exp(-risk))))
}

merge.cluster <- function(clusters, center, window.size) {
  copy.clusters <- clusters
  copy.clusters[which(copy.clusters %in% (center - c(-(window.size:1), 1:window.size)))] <- center
  return(sapply(copy.clusters, function(c) match(c, sort(unique(copy.clusters)))))
}

select.repr <- function(pvals, cst) {
  ind.repr <- sapply(1:max(cst), function(c) {
    cluster_elements <- cst == c
    top_within <- which.min(pvals[cluster_elements])
    if (length(top_within) == 0) top_within <- 1
    which(cluster_elements)[top_within]
  })

  clusters <- order(ind.repr)[cst]
  names(clusters) <- names(pvals.screen)
  ind.repr <- sort(ind.repr)

  return(list(clusters = clusters, ind.repr = ind.repr))
}
