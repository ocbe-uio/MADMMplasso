quick_func <- function(xz = NULL, xn) {
  as.vector(xz[1:xn] %o% xz[-(1:xn)])
}
