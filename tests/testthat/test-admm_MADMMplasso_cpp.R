# Auxiliary funcions ===========================================================
model <- function(beta0, theta0, beta, theta, X, Z) {
  N <- nrow(X)
  intercepts <- (matrix(1, N)) %*% beta0 + Z %*% ((theta0))
  shared_model <- X %*% (beta)
  intercepts + shared_model
}

reg_temp <- function(r, Z) {
  K <- ncol(Z)
  beta01 <- matrix(0, 1, ncol(r))
  theta01 <- matrix(0, ncol(Z), ncol(r))
  for (e in seq_len(ncol(r))) {
    my_one <- matrix(1, nrow(Z))
    my_w <- data.frame(Z, my_one)
    my_w <- as.matrix(my_w)
    my_inv <- solve(t(my_w) %*% my_w) # replaced pracma::pinv() to eliminate dep
    my_res <- my_inv %*% (t(my_w) %*% r[, e])
    beta01[e] <- matrix(my_res[(K + 1)])
    theta01[, e] <- matrix(my_res[1:K])
  }
  list(beta0 = beta01, theta0 = theta01)
}

# Generate the data ============================================================
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
e <- matrix(1, N)
e <- matrix(1, N)
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
e <- mvrnorm(N, mu = rep(0, 6), Sigma = esd)
y_train <- X %*% Beta + pliable + e
y <- y_train
colnames(y) <- 1:6
colnames(y) <- c(paste0("y", seq_len(ncol(y))))
TT <- tree_parms(y)
C <- TT$Tree
CW <- TT$Tw
N <- nrow(X)
p <- ncol(X)
K <- ncol(Z)
D <- ncol(y)
lambda <- rep(0.5, 6)
alpha <- 0.5
e.abs <- 1E-4
e.rel <- 1E-2
alpha <- 0.2
tol <- 1E-3
alph <- 1
rho <- 5
gg <- c(0.02, 0.02)
max_it <- 5000
my_W_hat <- generate_my_w(X = X, Z = Z)
svd.w <- svd(my_W_hat)
svd.w$tu <- t(svd.w$u)
svd.w$tv <- t(svd.w$v)
input <- 1:(ncol(y) * nrow(C))
multiple_of_D <- (input %% ncol(y)) == 0
I <- matrix(0, nrow = nrow(C) * ncol(y), ncol = ncol(y))
II <- input[multiple_of_D]
diag(I[seq_len(ncol(y)), ]) <- C[1, ] * (CW[1])
c_count <- 2
for (e in II[-length(II)]) {
  diag(I[c((e + 1):(c_count * ncol(y))), ]) <- C[c_count, ] * (CW[c_count])
  c_count <- 1 + c_count
}
new_I <- diag(t(I) %*% I)
new_G <- matrix(0, (p + p * K))
new_G[1:p] <- 1
new_G[-1:-p] <- 2
new_G[1:p] <- rho * (1 + new_G[1:p])
new_G[-1:-p] <- rho * (1 + new_G[-1:-p])
invmat <- list() # denominator of the beta estimates
for (rr in 1:D) {
  DD1 <- rho * (new_I[rr] + 1)
  DD2 <- new_G + DD1
  invmat[[rr]] <- DD2
}
beta0 <- matrix(0, 1, D) # estimates$Beta0
theta0 <- matrix(0, K, D)
beta <- (matrix(0, p, D))
beta_hat <- (matrix(0, p + p * (K), D))
V <- (array(0, c(p, 2 * (1 + K), D)))
O <- (array(0, c(p, 2 * (1 + K), D)))
E <- (matrix(0, ncol(y) * nrow(C), (p + p * K))) # response auxiliary
EE <- (array(0, c(p, (1 + K), D)))

# auxiliary variables for the L1 norm####
theta <- (array(0, c(p, K, D)))
Q <- (array(0, c(p, (1 + K), D)))
P <- (array(0, c(p, (1 + K), D)))
H <- (matrix(0, ncol(y) * nrow(C), (p + p * K))) # response multiplier
HH <- (array(0, c(p, (1 + K), D)))
r_current <- y
b <- reg_temp(r_current, Z) # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
beta0 <- b$beta0
theta0 <- b$theta0
new_y <- y - (matrix(1, N) %*% beta0 + Z %*% ((theta0)))
XtY <- crossprod((my_W_hat), (new_y))

# Testing the R function =======================================================
mprt <- FALSE
my_values <- suppressWarnings(suppressMessages(admm_MADMMplasso(
  beta0 = beta0, theta0 = theta0, beta = beta, beta_hat = beta_hat,
  theta = theta, rho, X, Z, max_it, W_hat = my_W_hat, XtY, y, N, e.abs,
  e.rel, alpha, lambda = lambda, alph, svd.w = svd.w, tree = TT,
  my_print = mprt, invmat = invmat, gg = gg
)))
beta <- my_values$beta
theta <- my_values$theta
converge <- my_values$converge
beta0 <- my_values$beta0
theta0 <- my_values$theta0 ### iteration
beta_hat <- my_values$beta_hat
y_hat <- my_values$y_hat

test_that("final objects have correct dimensions", {
  expect_identical(length(beta0), 6L)
  expect_identical(dim(theta0), c(4L, 6L))
  expect_identical(dim(beta), c(50L, 6L))
  expect_identical(dim(theta), c(50L, 4L, 6L))
  expect_length(converge, 1L)
  expect_identical(dim(beta_hat), c(250L, 6L))
  expect_identical(dim(y_hat), c(100L, 6L))
})

test_that("mean values of final objects are expected", {
  tole <- 1e-1
  expect_equal(mean(beta0), 5.132656e-02, tolerance = tole)
  expect_equal(mean(theta0), 5.123034e-02, tolerance = tole)
  expect_equal(mean(beta), 2.104393e-02, tolerance = tole)
  expect_equal(mean(theta), 2.841666e-04, tolerance = tole)
  expect_true(converge)
  expect_equal(mean(beta_hat), 4.436118e-03, tolerance = tole)
  expect_equal(mean(y_hat), -8.380419e-02, tolerance = tole)
})

# Testing the C++ function =====================================================
my_values_cpp <- admm_MADMMplasso_cpp(
  beta0, theta0, beta, beta_hat, theta, rho, X, Z, max_it, my_W_hat, XtY, y, N,
  e.abs, e.rel, alpha, lambda, alph, t(svd.w$u), t(svd.w$v), svd.w$d, TT$Tree,
  TT$Tw, gg, mprt
)

test_that("C++ function output structure", {
  expect_identical(length(my_values_cpp), length(my_values) - 1L)
})

test_that("Values are the same", {
  tl <- 1e-1
  expect_equal(my_values$beta0, my_values_cpp[[1]][, 1, 1], tolerance = tl)
  expect_equal(my_values$theta0, my_values_cpp[[2]][, , 1], tolerance = tl)
  expect_equal(my_values$beta, my_values_cpp[[3]][, , 1], tolerance = tl)
  expect_equal(my_values$theta, my_values_cpp[[4]], tolerance = tl)
  expect_equal(my_values$beta_hat, my_values_cpp[[6]][, , 1], tolerance = tl)
  expect_equal(my_values$y_hat, my_values_cpp[[7]][, , 1], tolerance = tl)
})
