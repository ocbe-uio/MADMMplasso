#' Generate the matrix W as used in Appendix I for use in the function.
#' @param X TODO: fill in description
#' @param Z TODO: fill in description
#' @param quad TODO: fill in description
#' @export
generate_my_w <- function(X = matrix(), Z = matrix(), quad = TRUE) {
  p1 <- ncol(X)
  p2 <- ncol(Z)

  if (quad == FALSE) {
    p <- p1
    if (p1 != p2) stop("To remove quadtratic terms p1 must be equal to p2")
    ind <- (1:p) * p + (2 * (1:p) + 1)
  }

  # Just in case we have only one oberservation? Not sure why I did this
  if (is.vector(X)) p1 <- length(X)
  if (is.vector(Z)) p2 <- length(Z)

  # Add the intercept
  x <- X
  z <- cbind(1, Z)

  W <- t(apply(cbind(x, z), 1, quick_func, xn = p1))
  if (quad == FALSE) {
    W <- W[, -ind]
  }

  return(W)
}
