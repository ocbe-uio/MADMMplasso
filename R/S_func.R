S_func <- function(x, a) { # Soft Thresholding Operator
  return(pmax(abs(x) - a, 0) * sign(x))
}
