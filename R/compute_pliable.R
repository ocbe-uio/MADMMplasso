#' @title Compute the interaction part of the model.
#' @description  Compute the interaction part of the model.
#' @param X N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical variables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param theta    theta coefficients for a single response ncol(X) by ncol(Z)
#' @return a vector of length N of the calculated interaction term for a single response
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

  xz_term
}
