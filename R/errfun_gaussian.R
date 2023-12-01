errfun_gaussian <- function(y, yhat, w = rep(1, length(y))) {
  (w * (y - yhat)^2)
}
