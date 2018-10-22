forward.sample <- function(x, pInit, Q, pEmit) {
  matOne <- array(1, dim = rep(length(pInit), 2))
  pObs <- log(pInit) + log(pEmit[1, x[1] + 1, ])
  tmp <- rep(0, length(pObs))

  for (i in 2:length(x)) {
    tmp <- matrixStats::rowLogSumExps(log(t(Q[i - 1, , ])) + matrix(rep(pObs, each = length(pInit)), nrow = length(pInit))) +
      log(pEmit[i, x[i] + 1, ])
    pObs <- tmp
  }
  return(matrixStats::logSumExp(pObs))
}
