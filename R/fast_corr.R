fast_corr <- function(A) {
  C <- crossprod(scale(A)) / (nrow(A) - 1)
  return(C)
}
