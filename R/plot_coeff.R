plot_coeff <- function(beta, theta, error, nz, p, K, D, nlambda, Lambda) {
  gg <- nlambda
  my_beta1 <- array(NA, c(p, nlambda, D))
  for (r in 1:nlambda) {
    for (i in 1:D) {
      my_beta1[, r, i] <- beta[[r]][, i]
    }
  }

  for (ii in 1:D) {
    my_beta <- my_beta1[, , ii]

    b <- apply(abs(my_beta), 2, sum)
    b <- log(unlist(Lambda[, ii]))
    n <- dim(my_beta)[2]
    matplot(b, t(my_beta), type = "n", col = "red", ylim = range(my_beta), xlab = "Log Lambda", ylab = (paste("coefficient", ii)))
    axis(
      side = 3, at = (as.matrix(b)), labels = paste(as.matrix(nz[1:gg])),
      tick = FALSE, line = 0
    )

    for (i in 1:(p)) {
      lines(b, (my_beta[i, ]), col = i + 1, lty = 1)

      text((min(b - 0.1)), my_beta[i, n], labels = i, cex = 0.7)
    }

    my_beta <- (my_beta)
    act <- which(rowSums(abs(my_beta)) > 0)
    theta1 <- array(0, c(p, K, (nlambda)))
    for (i in 1:(nlambda)) {
      theta1[, , i] <- matrix(unlist(theta[[i]][, , ii]), p, K)
    }

    ntheta <- apply(abs(theta1) > 1E-3, c(1, 3), sum)
    index <- b
    sbeta <- (my_beta)
    for (j in act) {
      for (i in seq_along(index)) {
        if (ntheta[j, i] > 0) text(index[i], sbeta[j, i], label = "x", cex = 0.7)
      }
    }
  }
}
