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
