# Auxiliary funcions ===========================================================
# TODO: ask why these are not in the package?
model <- function(beta0, theta0, beta, theta, X, Z) {
  p <- ncol(X)
  N <- nrow(X)
  K <- ncol(Z)
  D <- dim(beta0)[2]
  intercepts <- (matrix(1, N)) %*% beta0 + Z %*% ((theta0))
  shared_model <- X %*% (beta)
  pliable <- matrix(0, N, D)
  return(intercepts + shared_model)
}

reg <- function(r, Z) {
  K <- ncol(Z)
  beta01 <- matrix(0, 1, ncol(r))
  theta01 <- matrix(0, ncol(Z), ncol(r))
  for (e in 1:ncol(r)) {
    my_one <- matrix(1, nrow(Z))
    my_w <- data.frame(Z, my_one)
    my_w <- as.matrix(my_w)
    my_inv <- pracma::pinv(t(my_w) %*% my_w)
    my_res <- my_inv %*% (t(my_w) %*% r[, e])
    beta01[e] <- matrix(my_res[(K + 1)])
    theta01[, e] <- matrix(my_res[c(1:(K))])
  }
  return(list(beta0 = beta01, theta0 = theta01))
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
colnames(Beta) <- c(1:6)
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
colnames(y) <- c(1:6)
colnames(y) <- c(paste("y", 1:(ncol(y)), sep = ""))
TT <- tree.parms(y)
C <- TT$Tree
CW <- TT$Tw
N <- nrow(X)
p <- ncol(X)
K <- ncol(Z)
D <- dim(y)[2]
lambda <- rep(0.5, 6)
alpha <- 0.5
e.abs <- 1E-4
e.rel <- 1E-2
alpha <- .2
tol <- 1E-3
alph <- 1
rho <- 5
gg <- c(0.02, 0.02)
max_it <- 5000
my_W_hat <- generate.my.w(X = X, Z = Z, quad = TRUE)
svd.w <- svd(my_W_hat)
svd.w$tu <- t(svd.w$u)
svd.w$tv <- t(svd.w$v)
input <- 1:(dim(y)[2] * nrow(C))
multiple_of_D <- (input %% dim(y)[2]) == 0
I <- matrix(0, nrow = nrow(C) * dim(y)[2], ncol = dim(y)[2])
II <- input[multiple_of_D]
diag(I[c(1:dim(y)[2]), ]) <- C[1, ] * (CW[1])
c_count <- 2
for (e in II[-length(II)]) {
  diag(I[c((e + 1):(c_count * dim(y)[2])), ]) <- C[c_count, ] * (CW[c_count])
  c_count <- 1 + c_count
}
new_I <- diag(t(I) %*% I)
new_G <- matrix(0, (p + p * K))
new_G[c(1:p)] <- 1
new_G[-c(1:p)] <- 2
new_G[c(1:p)] <- rho * (1 + new_G[c(1:p)])
new_G[-c(1:p)] <- rho * (1 + new_G[-c(1:p)])
invmat <- list() # denominator of the beta estimates
for (rr in 1:D) {
  DD1 <- rho * (new_I[rr] + 1)
  DD2 <- new_G + DD1
  invmat[[rr]] <- DD2 # Matrix::chol2inv( Matrix::chol(new_sparse) )
}
beta0 <- matrix(0, 1, D) # estimates$Beta0
theta0 <- matrix(0, K, D)
beta <- (matrix(0, p, D))
beta_hat <- (matrix(0, p + p * (K), D))
V <- (array(0, c(p, 2 * (1 + K), D)))
O <- (array(0, c(p, 2 * (1 + K), D)))
E <- (matrix(0, dim(y)[2] * nrow(C), (p + p * K))) # response auxiliary
EE <- (array(0, c(p, (1 + K), D)))

# auxiliary variables for the L1 norm####
theta <- (array(0, c(p, K, D)))
Q <- (array(0, c(p, (1 + K), D)))
P <- (array(0, c(p, (1 + K), D)))
H <- (matrix(0, dim(y)[2] * nrow(C), (p + p * K))) # response multiplier
HH <- (array(0, c(p, (1 + K), D)))
r_current <- y #-model(beta0,theta0,beta=beta_hat, theta, X=W_hat, Z)
b <- reg(r_current, Z) # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
beta0 <- b$beta0
theta0 <- b$theta0
new_y <- y - (matrix(1, N) %*% beta0 + Z %*% ((theta0)))
XtY <- crossprod((my_W_hat), (new_y))
my_values <- suppressMessages(admm.MADMMplasso(
  beta0 = beta0, theta0 = theta0, beta = beta, beta_hat = beta_hat,
  theta = theta, rho, X, Z, max_it, W_hat = my_W_hat, XtY, y, N, p, K, e.abs,
  e.rel, alpha, lambda = lambda, alph, svd.w = svd.w, tree = TT,
  my_print = FALSE, invmat = invmat, V = V, Q = Q, E = E, EE = EE, O = O, P = P,
  H = H, HH = HH, cv = cv, gg = gg
))
beta <- my_values$beta
theta <- my_values$theta
converge <- my_values$converge
beta0 <- my_values$beta0
theta0 <- my_values$theta0 ### iteration
V <- my_values$V
Q <- my_values$Q
O <- my_values$O
P <- my_values$P
E <- my_values$E
H <- my_values$H
beta_hat <- my_values$beta_hat
y_hat <- my_values$y_hat

test_that("final objects have correct dimensions", {
	expect_equal(dim(beta0), c(1, 6))
	expect_equal(dim(theta0), c(4, 6))
	expect_equal(dim(beta), c(50, 6))
	expect_equal(dim(theta), c(50, 4, 6))
	expect_equal(length(converge), 1)
	expect_equal(dim(V), c(50, 10, 6))
	expect_equal(dim(Q), c(50, 5, 6))
	expect_equal(dim(O), c(50, 10, 6))
	expect_equal(dim(P), c(50, 5, 6))
	expect_equal(dim(E), c(18, 250))
	expect_equal(dim(H), c(18, 250))
	expect_equal(dim(EE), c(50, 5, 6))
	expect_equal(dim(HH), c(50, 5, 6))
	expect_equal(dim(beta_hat), c(250, 6))
	expect_equal(dim(y_hat), c(100, 6))
})

test_that("mean values of final objects are expected", {
  tolerance <- 1e-4
	expect_equal(mean(beta0), 5.132656e-02, tolerance = tolerance)
	expect_equal(mean(theta0), 5.123034e-02, tolerance = tolerance)
	expect_equal(mean(beta), 2.104393e-02, tolerance = tolerance)
	expect_equal(mean(theta), 2.841666e-04, tolerance = tolerance)
	expect_equal(converge, TRUE)
	expect_equal(mean(V), 2.329982e-03, tolerance = tolerance)
	expect_equal(mean(Q), 4.439442e-03, tolerance = tolerance)
	expect_equal(mean(O), 1.462051e-03, tolerance = tolerance)
	expect_equal(mean(P), 6.868194e-04, tolerance = tolerance)
	expect_equal(mean(E), 1.477602e-03, tolerance = tolerance)
	expect_equal(mean(H), -2.318846e-05, tolerance = tolerance)
	expect_equal(mean(EE), 0, tolerance = tolerance)
	expect_equal(mean(HH), -6.948488e-05, tolerance = tolerance)
	expect_equal(mean(beta_hat), 4.436118e-03, tolerance = tolerance)
	expect_equal(mean(y_hat), -8.380419e-02, tolerance = tolerance)
})

# Testing the C++ function =====================================================
my_values_cpp <- admm_MADMMplasso_cpp(
  beta0 = beta0, theta0 = theta0, beta = beta, beta_hat = beta_hat,
  theta = theta, rho, X, Z, max_it, W_hat = my_W_hat, XtY, y, N, p, K, e.abs,
  e.rel, alpha, lambda = lambda, alph, svd_w = svd.w, tree = TT,
  my_print = FALSE, invmat = invmat, V = V, Q = Q, E = E, EE = EE, O = O, P = P,
  H = H, HH = HH, gg = gg
)

test_that("C++ function output structure", {
  expect_equal(length(my_values_cpp), 16)
})
