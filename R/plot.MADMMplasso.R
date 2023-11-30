#' @export
plot.MADMMplasso <- function(x, ...) {
  fit <- x
  beta <- fit$beta
  theta <- fit$theta
  p <- nrow(fit$beta[[1]])
  K <- nrow(fit$theta0[[1]])
  D <- ncol(fit$theta0[[1]])
  nlambda <- length(fit$Lambdas[, 1])

  plot_coeff(beta, theta, error = fit$path$OBJ_main, nz = fit$non_zero, p, K, D, nlambda, Lambda = fit$Lambdas)
}
