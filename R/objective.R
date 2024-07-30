objective <- function(beta0, theta0, beta, theta, X, Z, y, alpha, lambda, p, N, IB, W, beta1) {
  loss <- (norm(y - model_p(beta0, theta0, beta = beta1, X = W, Z), type = "F")^2)
  mse <- (1 / (2 * N)) * loss

  l_1 <- sum(abs(beta))
  pliable_norm <- matrix(0, dim(y)[2])
  for (ee in 1:dim(y)[2]) {
    beta11 <- beta[, ee]
    theta11 <- theta[, , ee]
    norm_1_l <- lapply(
      seq_len(p),
      function(g) (lambda[ee] * (1 - alpha) * (norm(matrix(c(beta11[g], theta11[g, ])), type = "F") + norm(matrix(c(theta11[g, ])), type = "F")) + lambda[ee] * alpha * sum(abs(theta11[g, ])))
    )

    pliable_norm[ee] <- sum(unlist(norm_1_l))
  }

  # Output
  mse + (1 - alpha) * min(lambda / 4) * IB + (1 - alpha) * min(lambda / 4) * l_1 + sum(pliable_norm)
}
