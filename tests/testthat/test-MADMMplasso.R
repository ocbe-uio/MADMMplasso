# Setting up objects =========================================================
set.seed(1235)
N <- 100
p <- 50
nz <- 4
K <- nz
X <- matrix(rnorm(n = N * p), nrow = N, ncol = p)
mx <- colMeans(X)
sx <- sqrt(apply(X, 2, var))
X <- scale(X, mx, sx)
X <- matrix(as.numeric(X), N, p)
Z <- matrix(rnorm(N * nz), N, nz)
mz <- colMeans(Z)
sz <- sqrt(apply(Z, 2, var))
Z <- scale(Z, mz, sz)
beta_1 <- rep(x = 0, times = p)
beta_2 <- rep(x = 0, times = p)
beta_3 <- rep(x = 0, times = p)
beta_4 <- rep(x = 0, times = p)
beta_5 <- rep(x = 0, times = p)
beta_6 <- rep(x = 0, times = p)

beta_1[1:5] <- c(2, 2, 2, 2, 2)
beta_2[1:5] <- c(2, 2, 2, 2, 2)
beta_3[6:10] <- c(2, 2, 2, -2, -2)
beta_4[6:10] <- c(2, 2, 2, -2, -2)
beta_5[11:15] <- c(-2, -2, -2, -2, -2)
beta_6[11:15] <- c(-2, -2, -2, -2, -2)

Beta <- cbind(beta_1, beta_2, beta_3, beta_4, beta_5, beta_6)
colnames(Beta) <- 1:6

theta <- array(0, c(p, K, 6))
theta[1, 1, 1] <- 2
theta[3, 2, 1] <- 2
theta[4, 3, 1] <- -2
theta[5, 4, 1] <- -2
theta[1, 1, 2] <- 2
theta[3, 2, 2] <- 2
theta[4, 3, 2] <- -2
theta[5, 4, 2] <- -2
theta[6, 1, 3] <- 2
theta[8, 2, 3] <- 2
theta[9, 3, 3] <- -2
theta[10, 4, 3] <- -2
theta[6, 1, 4] <- 2
theta[8, 2, 4] <- 2
theta[9, 3, 4] <- -2
theta[10, 4, 4] <- -2
theta[11, 1, 5] <- 2
theta[13, 2, 5] <- 2
theta[14, 3, 5] <- -2
theta[15, 4, 5] <- -2
theta[11, 1, 6] <- 2
theta[13, 2, 6] <- 2
theta[14, 3, 6] <- -2
theta[15, 4, 6] <- -2

pliable <- matrix(0, N, 6)
for (e in 1:6) {
  pliable[, e] <- compute_pliable(X, Z, theta[, , e])
}

esd <- diag(6)
e <- MASS::mvrnorm(N, mu = rep(0, 6), Sigma = esd)
y_train <- X %*% Beta + pliable + e
y <- y_train

colnames(y) <- c(paste0("y", seq_len(ncol(y))))
TT <- tree_parms(y)
gg1 <- matrix(0, 2, 2)
gg1[1, ] <- c(0.02, 0.02)
gg1[2, ] <- c(0.02, 0.02)

nlambda <- 1
e.abs <- 1E-4
e.rel <- 1E-2
alpha <- 0.2
tol <- 1E-3

# Running MADMMplasso ========================================================
set.seed(9356219)
fit_C <- MADMMplasso(
  X, Z, y,
  alpha = alpha, my_lambda = matrix(rep(0.2, ncol(y)), 1),
  lambda_min = 0.001, max_it = 5000, e.abs = e.abs, e.rel = e.rel, maxgrid = nlambda,
  nlambda = nlambda, rho = 5, tree = TT, my_print = FALSE, alph = 1,
  pal = TRUE, gg = gg1, tol = tol, cl = 1
)
set.seed(9356219)
fit_R <- suppressWarnings(
  suppressMessages(
    MADMMplasso(
      X, Z, y,
      alpha = alpha, my_lambda = matrix(rep(0.2, ncol(y)), 1),
      lambda_min = 0.001, max_it = 5000, e.abs = e.abs, e.rel = e.rel, maxgrid = nlambda,
      nlambda = nlambda, rho = 5, tree = TT, my_print = FALSE, alph = 1,
      pal = TRUE, gg = gg1, tol = tol, cl = 1, legacy = TRUE
    )
  )
)

test_that("C++ and R versions basically output the same thing", {
  expect_named(fit_C$beta, names(fit_R$beta))
  tl <- 1e1
  expect_equal(fit_C$beta0[[1]], as.matrix(fit_R$beta0[[1]]), tolerance = tl)
  expect_equal(as.vector(fit_C$beta[[1]]), as.vector(fit_R$beta[[1]]), tolerance = tl)
  expect_equal(as.vector(fit_C$BETA_hat[[1]]), as.vector(fit_R$BETA_hat[[1]]), tolerance = tl)
  expect_equal(fit_C$theta0[[1]], fit_R$theta0[[1]], tolerance = tl)
  for (i in 1:6) {
    expect_equal(
      as.vector(fit_C$theta[[1]][, , i]),
      as.vector(fit_R$theta[[1]][, , i]),
      tolerance = tl
    )
  }
  expect_equal(fit_C$path, fit_R$path, tolerance = tl)
  expect_identical(fit_C$Lambdas, fit_R$Lambdas)
  expect_identical(fit_C$non_zero[1], fit_R$non_zero)
  expect_identical(fit_C$LOSS[1], fit_R$LOSS)
  expect_equal(fit_C$Y_HAT[[1]], fit_R$Y_HAT[[1]], tolerance = tl)
  expect_identical(fit_C$gg, fit_R$gg)
})
