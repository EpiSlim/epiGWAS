forward.sample <- function(x, pInit, Q, pEmit) {
  matOne <- array(1, dim = rep(length(pInit), 2))
  pObs <- log(pInit) + log(pEmit[1, x[1] + 1, ])
  for (i in 2:length(x)) {
    pObs <- matrixStats::rowLogSumExps(log(t(Q[i - 1, , ])) + matOne %*% diag(pObs)) +
      log(pEmit[i, x[i] + 1, ])
  }
  return(matrixStats::logSumExp(pObs))
}

forward <- function(X, pInit, Q, pEmit) {
  # Modify the input controls
  stopifnot(length(pInit) == dim(Q)[2])
  stopifnot(dim(pEmit)[3] == dim(Q)[2])
  stopifnot(dim(pEmit)[1] == dim(Q)[1] + 1)
  stopifnot(dim(Q)[2] == dim(Q)[3])

  `%operator%` <- foreach::`%dopar%`

  pObs <- foreach::foreach(i = 1:nrow(X), .export = c("forward.sample"), .combine = "c") %operator% {
    return(forward.sample(X[i, ], pInit, Q, pEmit))
  }

  return(pObs)
}

cond.prob <- function(X, target.name, out.path=NULL, fp.path="bin/fastPHASE", n.state=12, n.iter=25) {
  # Running fastPHASE
  Xinp_file <- SNPknock::SNPknock.fp.writeX(X)
  fp.outPath <- SNPknock::SNPknock.fp.runFastPhase(fp.path, Xinp_file, K = n.state, numit = n.iter, out_path = out.path)

  # Loading the fitted the Hidden Markov Model
  r_file <- paste(fp.outPath, "_rhat.txt", sep = "")
  theta_file <- paste(fp.outPath, "_thetahat.txt", sep = "")
  alpha_file <- paste(fp.outPath, "_alphahat.txt", sep = "")
  char_file <- paste(fp.outPath, "_origchars", sep = "")

  hmm <- SNPknock::SNPknock.fp.loadFit(r_file, theta_file, alpha_file, X[1, ])

  # Computing the propensity scores with the forward algorithm
  pObs <- matrix(nrow = nrow(X), ncol = 3)
  copy.X <- X
  copy.X[, target.name] <- 0
  pObs[, 1] <- forward(copy.X, hmm$pInit, hmm$Q, hmm$pEmit)
  copy.X[, target.name] <- 1
  pObs[, 2] <- forward(copy.X, hmm$pInit, hmm$Q, hmm$pEmit)
  copy.X[, target.name] <- 2
  pObs[, 3] <- forward(copy.X, hmm$pInit, hmm$Q, hmm$pEmit)
  pObs <- t(apply(pObs, 1, function(x) return(x - matrixStats::logSumExp(x))))

  return(cbind(exp(pObs[, 1]), exp(pObs[, 2]) + exp(pObs[, 3])))
}
