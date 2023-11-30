#' @title TODO: fill this field
#' @description TODO: fill this field
#' @param X TODO: fill this field
#' @param Z TODO: fill this field
#' @param theta TODO: fill this field
#' @return TODO: fill this field
#' @export
compute_pliable <- function(X, Z, theta) {
  p <- ncol(X)
  N <- nrow(X)
  K <- ncol(Z)

  xz_theta <- lapply(
    seq_len(p),
    function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j]
  )
  xz_term <- (Reduce(f = "+", x = xz_theta))

  return(xz_term)
}
