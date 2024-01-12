#' @title Simulate data for the model. This is the second simulation data used in the paper
#' @description Simulate data for the model
#' @param p column for X which is the main effect
#' @param n number of observations
#' @param m number of responses
#' @param nz  number of modifiers
#' @param rho  values to be used in the covariance matrix when generating  X
#' @param B.elem  the value of the non-zero elements in beta?
#' @return  The simulated data with the following components:
#'  Beta: matrix of actual beta coefficients  p by m
#'  Theta: a m by p by K array of actual theta coefficients
#'  Y: a N by m matrix of response variables
#'  X: a N by p matrix of covariates
#'  Z: a N by K matrix of modifiers

#' @export
sim2 <- function(p = 500, n = 100, m = 24, nz = 4, rho = 0.4, B.elem = 0.5) {
  b <- 10
  if (!is.na(p)) {
    # generate covariance matrix
    Ap1 <- matrix(rep(rho, (p / b)^2), nrow = p / b)
    diag(Ap1) <- rep(1, p / b)
    Xsigma1 <- Ap1
    for (i in 2:b) {
      Xsigma1 <- bdiag(Xsigma1, Ap1)
    }

    Xsigma <- rbind(cbind(Xsigma1))
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Xsigma)
    esd <- diag(m)
    e <- MASS::mvrnorm(n, mu = rep(0, m), Sigma = esd)

    ## generate beta1 matrix
    Beta1 <- matrix(0, nrow = m, ncol = p)
    theta <- array(0, c(p, nz, m))
    Beta1[, 1] <- B.elem

    for (i in 1:2) {
      Beta1[((i - 1) * m / 2 + 1):(i * m / 2), (1 + (i - 1) * 2 + 1):(1 + i * 2)] <- B.elem
    }
    for (i in 1:4) {
      Beta1[((i - 1) * m / 4 + 1):(i * m / 4), (1 + 4 * 2 + (i - 1) * 4 + 1):(1 + 4 * 2 + i * 4)] <- B.elem
    }
    for (i in 1:8) {
      Beta1[((i - 1) * m / 8 + 1):(i * m / 8), (1 + 2 * 2 + 4 * 4 + (i - 1) * 8 + 1):(1 + 2 * 2 + 4 * 4 + i * 8)] <- B.elem
    }

    theta[30, 1, 1] <- 0.6
    theta[31, 2, 1] <- 0.6
    theta[32, 3, 1] <- -0.6
    theta[33, 4, 1] <- -0.6
    theta[30, 1, 2] <- 0.6
    theta[31, 2, 2] <- 0.6
    theta[32, 3, 2] <- -0.6
    theta[33, 4, 2] <- -0.6

    theta[35, 1, 5] <- 0.6
    theta[36, 2, 5] <- 0.6
    theta[37, 3, 5] <- -0.6
    theta[38, 4, 5] <- -0.6
    theta[35, 1, 6] <- 0.6
    theta[36, 2, 6] <- 0.6
    theta[37, 3, 6] <- -0.6
    theta[38, 4, 6] <- -0.6

    theta[40, 1, 8] <- 0.6
    theta[41, 2, 8] <- 0.6
    theta[42, 3, 8] <- -0.6
    theta[43, 4, 8] <- -0.6
    theta[40, 1, 9] <- 0.6
    theta[41, 2, 9] <- 0.6
    theta[42, 3, 9] <- -0.6
    theta[43, 4, 9] <- -0.6

    theta[48, 1, 10] <- 0.6
    theta[49, 2, 10] <- 0.6
    theta[50, 3, 10] <- -0.6
    theta[51, 4, 10] <- -0.6
    theta[48, 1, 11] <- 0.6
    theta[49, 2, 11] <- 0.6
    theta[50, 3, 11] <- -0.6
    theta[51, 4, 11] <- -0.6

    theta[57, 1, 13] <- 0.6
    theta[58, 2, 15] <- 0.6
    theta[59, 3, 15] <- -0.6
    theta[60, 4, 15] <- -0.6
    theta[57, 1, 15] <- 0.6
    theta[58, 2, 17] <- 0.6
    theta[59, 3, 17] <- -0.6
    theta[60, 4, 17] <- -0.6

    theta[63, 1, 16] <- 0.6
    theta[64, 2, 16] <- 0.6
    theta[65, 3, 16] <- -0.6
    theta[66, 4, 16] <- -0.6
    theta[63, 1, 17] <- 0.6
    theta[64, 2, 17] <- 0.6
    theta[65, 3, 17] <- -0.6
    theta[66, 4, 17] <- -0.6

    theta[80, 1, 21] <- 0.6
    theta[81, 2, 21] <- 0.6
    theta[82, 3, 21] <- -0.6
    theta[83, 4, 21] <- -0.6
    theta[80, 1, 22] <- 0.6
    theta[81, 2, 22] <- 0.6
    theta[82, 3, 22] <- -0.6
    theta[83, 4, 22] <- -0.6

  }
  Beta <- t((Beta1))

  mx <- colMeans(X)

  sx <- sqrt(apply(X, 2, var))
  X <- scale(X, mx, sx)
  Z <- matrix(rbinom(n = n * nz, size = 1, prob = 0.5), nrow = n, ncol = nz)
  pliable <- matrix(0, n, m)
  for (ee in 1:m) {
    pliable[, ee] <- compute_pliable(X, Z, theta[, , ee])
  }
  Y <- X %*% Beta + pliable + e
  out <- list(Y = Y, X = X, Z = Z, Beta = Beta, Theta = theta, e = e, p = p)
  return(out)
}
