GBOOST <- function(X, A, Y) {
  `%operator%` <- foreach::`%dopar%`

  InteractionBOOST <- foreach::foreach(i = 1:ncol(X), .combine = "c") %operator% {
    cov.AX <- data.frame(x1 = factor(A), x2 = factor(X[, i]), y = Y)
    fit01 <- glm(y ~ x1 + x2, family = "binomial", data = cov.AX)
    fit2 <- glm(y ~ x1 + x2 + x1 * x2, family = "binomial", data = cov.AX)

    return(fit01$dev - fit2$dev)
  }

  return(InteractionBOOST)
}
