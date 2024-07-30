fast_corr <- function(A) {
  crossprod(scale(A)) / (nrow(A) - 1)
}
