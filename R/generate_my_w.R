
#' @title Generate the matrix W as seen in equation 8  for use in the function.
#' @description Generate the matrix W as seen in equation 8  for use in the function.
#' @param X N by p matrix of predictors
#' @param Z N by nz matrix of modifying variables. The elements of z
#' may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical variables should be coded by 0-1 dummy variables: for a k-level
#' variable, one can use either k or k-1  dummy variables.
#' @param quad TODO: fill in description
#' 
#' @return  Generated W matrix nrow(X) by (ncol(X)+ncol(X) by ncol(Z))
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
