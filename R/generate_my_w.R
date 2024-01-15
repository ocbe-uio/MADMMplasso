
#' @title Generate the matrix W as seen in equation 8  for use in the function.
#' @description Generate the matrix W as seen in equation 8  for use in the function.
#' @param X N by p matrix of predictors
#' @param Z N by nz matrix of modifying variables. The elements of z
#' may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical variables should be coded by 0-1 dummy variables: for a k-level
#' variable, one can use either k or k-1  dummy variables.
#'
#' @return  Generated W matrix nrow(X) by (ncol(X)+ncol(X) by ncol(Z))
#' @export
generate_my_w <- function(X = matrix(), Z = matrix()) {
  p1 <- ncol(X)

  # Just in case we have only one oberservation? Not sure why I did this
  if (is.vector(X)) p1 <- length(X)

  # Add the intercept
  x <- X
  z <- cbind(1, Z)

  W <- t(apply(cbind(x, z), 1, quick_func, xn = p1))

  return(W)
}
