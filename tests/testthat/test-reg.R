# Original function ============================================================
reg_R <- function(r, Z) {
  beta01 <- matrix(0, 1, ncol(r))
  theta01 <- matrix(0, ncol(Z), ncol(r))
  for (e in seq_len(ncol(r))) {
    new1 <- lm(r[, e] ~ Z, singular.ok = TRUE)
    beta01[e] <- matrix(new1$coefficients[1])
    theta01[, e] <- as.vector(new1$coefficients[-1])
  }
  list(beta0 = beta01, theta0 = theta01)
}

# Testing ======================================================================
reps <- 10L
n_obs <- rpois(reps, lambda = 10L)
n_vars <- sample(2:10, reps, replace = TRUE)
test_that("reg() produces the correct output", {
  for (rp in seq_len(reps)) {
    r <- matrix(rnorm(n_obs[rp] * n_vars[rp]), n_obs[rp], n_vars[rp])
    z <- as.matrix(sample(0:1, n_obs[rp], replace = TRUE))
    expect_identical(reg(r, z)[1, ], reg_R(r, z)[[1]][1, ], tolerance = 1e-10)
    expect_identical(reg(r, z)[2, ], reg_R(r, z)[[2]][1, ], tolerance = 1e-10)
  }
})
