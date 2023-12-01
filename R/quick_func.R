quick_func <- function(xz = c(), xn) {
  as.vector(xz[1:xn] %o% xz[-(1:xn)])
}
