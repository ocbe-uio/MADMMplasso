fast_corr <- function(A) {
  C <- crossprod(scale(A)) / (dim(A)[1] - 1)
  return(C)
}
