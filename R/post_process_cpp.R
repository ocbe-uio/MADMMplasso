post_process_cpp <- function(lst) {
  array2list <- function(ra) {
    apply(ra, 3, function(x) x, simplify = FALSE)
  }
  lst$BETA0 <- array2list(lst$BETA0)
  lst$THETA0 <- array2list(lst$THETA0)
  lst$BETA <- array2list(lst$BETA)
  lst$BETA_hat <- array2list(lst$BETA_hat)
  lst$Y_HAT <- array2list(lst$Y_HAT)
  lst
}
