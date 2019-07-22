num.infoIc <- function(coefficients, FUN, X, data,cublim)
{
  FUN <- get(FUN, inherits = TRUE)
  values <- FUN(coefficients, X, data,cublim)
  p <- length(values)
  Info <- matrix(0, p, p)
  h <- rep(0, p)
  delta <- cbind((abs(coefficients) + 1e-012) * 0.0001, rep(1e-012, p))
  delta <- apply(delta, 1, max)
  for(i in 1:p) {
    h[i] <- delta[i]
    new.values <- FUN(coefficients + h, X, data,cublim)
    Info[, i] <- (new.values - values)/delta[i]
    h[i] <- 0
  }
  Info
}