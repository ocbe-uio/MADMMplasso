count_nonzero_a <- function(x) {
  if (length(dim(x)) == 3) {
    count1 <- matrix(0, dim(x)[3])
    for (ww in seq_len(dim(x)[3])) {
      n <- sum(x[, , ww] != 0)
      count1[ww] <- n
    }
    n <- max(count1)
  } else {
    count1 <- matrix(0, ncol(x))
    for (ww in seq_len(ncol(x))) {
      n <- sum(x[, ww] != 0)
      count1[ww] <- n
    }
    n <- max(count1)
  }
  n
}
