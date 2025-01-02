#' @title Fit the ADMM part of  model for the given lambda values
#' @description This function fits a multi-response pliable lasso model over a path of regularization values.
#' @inheritParams MADMMplasso
#' @inheritParams cv_MADMMplasso
#' @param beta0 a vector of length ncol(y) of estimated beta_0 coefficients
#' @param theta0  matrix of the initial theta_0 coefficients  ncol(Z) by ncol(y)
#' @param beta  a matrix of the initial beta coefficients    ncol(X) by ncol(y)
#' @param beta_hat  a matrix of the initial beta and theta coefficients    (ncol(X)+ncol(X) by ncol(Z)) by ncol(y)
#' @param theta an array of initial theta coefficients    ncol(X) by ncol(Z) by ncol(y)
#' @param rho1 the Lagrange variable for the ADMM which is usually included as rho in the MADMMplasso call.
#' @param W_hat N by (p+(p by nz)) of the main and interaction predictors. This generated internally  when MADMMplasso is called or by using the function generate_my_w.
#' @param XtY a matrix formed by multiplying the transpose of X by y.
#' @param N nrow(X)
#' @param svd.w singular value decomposition of W
#' @param invmat A list of length ncol(y), each containing the C_d part of equation 32 in the paper
#' @param gg penalty terms for the tree structure for lambda_1 and  lambda_2 for the ADMM call.
#' @return  predicted values for the ADMM part

#' beta0:  estimated beta_0 coefficients  having a size of 1 by ncol(y)
#'
#' beta: estimated beta coefficients  having a matrix   ncol(X) by ncol(y)
#'
#' BETA_hat:  estimated beta and theta coefficients  having a matrix   (ncol(X)+ncol(X) by ncol(Z)) by ncol(y)
#'
#'  theta0:  estimated theta_0 coefficients  having a matrix   ncol(Z) by ncol(y)
#'
#'  theta:  estimated theta coefficients  having a an array   ncol(X) by ncol(Z) by ncol(y)
#'  converge: did the algorithm converge?
#'
#'  Y_HAT:  predicted response nrow(X) by ncol(y)



#' @export
admm_MADMMplasso <- function(beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, W_hat, XtY, y, N, e.abs, e.rel, alpha, lambda, alph, svd.w, tree, my_print, invmat, gg = 0.2) {
  TT <- tree

  C <- TT$Tree
  CW <- TT$Tw
  svd.w$tu <- t(svd.w$u)
  svd.w$tv <- t(svd.w$v)
  D <- ncol(y)
  p <- ncol(X)
  K <- ncol(Z)

  V <- (array(0, c(p, 2 * (1 + K), D)))
  O <- (array(0, c(p, 2 * (1 + K), D)))
  E <- (matrix(0, ncol(y) * nrow(C), (p + p * K))) # response auxiliary
  EE <- (array(0, c(p, (1 + K), D)))

  Q <- (array(0, c(p, (1 + K), D)))
  P <- (array(0, c(p, (1 + K), D)))
  H <- (matrix(0, ncol(y) * nrow(C), (p + p * K))) # response multiplier
  HH <- (array(0, c(p, (1 + K), D)))

  ### for response groups ###############################################################

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

  ###### overlapping group fro covariates########

  G <- matrix(0, 2 * (1 + K), (1 + K))
  diag_G <- matrix(0, (K + 1), (K + 1))
  diag(diag_G) <- 1
  for (i in 1:(K + 1)) {
    G[i, ] <- diag_G[i, ]
  }
  for (i in (K + 3):(2 * (K + 1))) {
    G[i, ] <- diag_G[(i - (K + 1)), ]
  }

  ####################################################################

  # auxiliary variables for the L1 norm####

  ################################################

  V_old <- V
  Q_old <- Q
  E_old <- E
  EE_old <- EE
  res_pri <- 0
  res_dual <- 0
  obj <- NULL

  SVD_D <- Diagonal(x = svd.w$d)
  R_svd <- (svd.w$u %*% SVD_D) / N
  R_svd_inv <- solve(R_svd)

  rho <- rho1
  Big_beta11 <- V
  for (i in 2:max_it) {
    r_current <- (y - model_intercept(beta_hat, W_hat))
    b <- reg(r_current, Z) # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
    beta0 <- b[1, ]
    theta0 <- b[-1, ]

    new_y <- y - (matrix(1, N) %*% beta0 + Z %*% ((theta0)))

    XtY <- crossprod((W_hat), (new_y))

    res_val <- rho * (t(I) %*% (E) - (t(I) %*% (H)))

    v.diff1 <- matrix(0, D)
    q.diff1 <- matrix(0, D)
    ee.diff1 <- matrix(0, D)

    new_G <- matrix(0, (p + p * K))
    new_G[1:p] <- 1
    new_G[-1:-p] <- 2
    new_G[1:p] <- rho * (1 + new_G[1:p])
    new_G[-1:-p] <- rho * (1 + new_G[-1:-p])

    invmat <- lapply(seq_len(D), function(j) new_G + rho * (new_I[j] + 1))

    for (jj in 1:D) {
      group <- (rho) * (t(G) %*% t(V[, , jj]) - t(G) %*% t(O[, , jj]))
      group1 <- group[1, ]
      group2 <- t(group[-1, ])
      new_group <- matrix(0, p, (K + 1))

      new_group[, 1] <- group1
      new_group[, -1] <- group2

      my_beta_jj <- XtY[, jj] / N + as.vector(new_group) + as.vector(res_val[jj, ]) + as.vector(rho * (Q[, , jj] - P[, , jj])) + as.vector(rho * (EE[, , jj] - HH[, , jj]))

      my_beta_jj <- matrix(my_beta_jj, ncol = 1)

      DD3 <- Diagonal(x = 1 / invmat[[jj]])

      part_z <- DD3 %*% t(W_hat)
      part_y <- DD3 %*% my_beta_jj

      beta_hat_j <- solve(R_svd_inv + svd.w$tv %*% part_z)
      beta_hat_j <- beta_hat_j %*% (svd.w$tv %*% part_y)
      beta_hat_j <- part_z %*% beta_hat_j
      beta_hat_JJ <- part_y - beta_hat_j
      beta_hat_JJ <- matrix(beta_hat_JJ, ncol = 1)

      beta_hat[, jj] <- beta_hat_JJ

      beta_hat1 <- matrix(beta_hat_JJ, p, (1 + K))

      b_hat <- alph * beta_hat1 + (1 - alph) * Q[, , jj]
      Q[, 1, jj] <- b_hat[, 1] + (P[, 1, jj])
      new.mat <- b_hat[, -1] + P[, -1, jj]
      Q[, -1, jj] <- sign(new.mat) * pmax(abs(new.mat) - ((alpha * lambda[jj]) / (rho)), 0)
      b_hat <- alph * beta_hat1 + (1 - alph) * EE[, , jj]
      new.mat <- b_hat + HH[, , jj]
      row.norm1 <- sqrt(apply(new.mat^2, 1, sum, na.rm = TRUE))
      coef.term1 <- pmax(1 - (gg[2]) / rho / (row.norm1), 0)
      ee1 <- scale(t(as.matrix(new.mat)), center = FALSE, scale = 1 / coef.term1)
      EE[, , jj] <- t(ee1)

      Big_beta <- t(tcrossprod(G, (beta_hat1)))
      Big_beta11[, , jj] <- Big_beta
      Big_beta1 <- alph * Big_beta + (1 - alph) * V[, , jj]

      # Now we have the main part.
      new.mat <- Big_beta1 + O[, , jj]
      new.mat1 <- new.mat[, 1:(K + 1)]
      new.mat2 <- new.mat[, -1:-(K + 1)]
      row.norm1 <- sqrt(apply(new.mat1^2, 1, sum, na.rm = TRUE))
      row.norm2 <- sqrt(apply(new.mat2^2, 1, sum, na.rm = TRUE))

      coef.term1 <- pmax(1 - ((1 - alpha) * lambda[jj]) / (rho) / (row.norm1), 0)
      coef.term2 <- pmax(1 - ((1 - alpha) * lambda[jj]) / (rho) / (row.norm2), 0)
      N_V1 <- scale(t(new.mat1), center = FALSE, scale = 1 / coef.term1)
      N_V2 <- scale(t(new.mat2), center = FALSE, scale = 1 / coef.term2)

      V[, , jj] <- cbind(t(N_V1), t(N_V2))

      P[, , jj] <- P[, , jj] + beta_hat1 - Q[, , jj]
      HH[, , jj] <- HH[, , jj] + beta_hat1 - EE[, , jj]
      O[, , jj] <- O[, , jj] + Big_beta - V[, , jj]

      v.diff1[jj] <- sum(((Big_beta - V[, , jj]))^2, na.rm = TRUE)
      q.diff1[jj] <- sum(((beta_hat1 - Q[, , jj]))^2, na.rm = TRUE)
      ee.diff1[jj] <- sum(((beta_hat1 - EE[, , jj]))^2, na.rm = TRUE)
    }

    ############ to estimate E  ##################
    Big_beta_respone <- ((I) %*% t(beta_hat))
    b_hat_response <- alph * Big_beta_respone + (1 - alph) * E
    new.mat <- b_hat_response + H

    new.mat_group <- (array(NA, c(p + p * K, ncol(y), nrow(C))))
    beta.group <- (array(NA, c(p + p * K, ncol(y), nrow(C))))
    N_E <- list()
    II <- input[multiple_of_D]
    new.mat_group[, , 1] <- t((new.mat[seq_len(ncol(y)), ]))
    beta.group[, , 1] <- t((Big_beta_respone[seq_len(ncol(y)), ]))

    beta_transform <- matrix(0, p, (K + 1) * ncol(y))
    beta_transform[, 1:(1 + K)] <- matrix(new.mat_group[, 1, 1], ncol = (K + 1), nrow = p)
    input2 <- 1:(ncol(y) * (1 + K))
    multiple_of_K <- (input2 %% (K + 1)) == 0
    II2 <- input2[multiple_of_K]
    e2 <- II2[-length(II2)][1]

    for (c_count2 in 2:ncol(y)) {
      beta_transform[, c((e2 + 1):(c_count2 * (1 + K)))] <- matrix(new.mat_group[, c_count2, 1], ncol = (K + 1), nrow = p)
      e2 <- II2[c_count2]
    }

    norm_res <- ((apply(beta_transform, 1, twonorm)))
    coef.term1 <- pmax(1 - (gg[1]) / rho / (norm_res), 0)

    N_E1 <- scale(t(beta_transform), center = FALSE, scale = 1 / coef.term1)

    N_E1 <- t(N_E1)
    beta_transform1 <- matrix(0, p + p * K, ncol(y))
    beta_transform1[, 1] <- as.vector(N_E1[, 1:(K + 1)])

    input3 <- 1:(ncol(y) * (1 + K))
    multiple_of_K <- (input3 %% (K + 1)) == 0
    II3 <- input3[multiple_of_K]
    e3 <- II3[-length(II3)][1]

    for (c_count3 in 2:ncol(y)) {
      beta_transform1[, c_count3] <- as.vector(N_E1[, c((e3 + 1):((K + 1) * c_count3))])
      e3 <- II3[c_count3]
    }

    N_E[[1]] <- (t(beta_transform1))

    e <- II[-length(II)][1]
    for (c_count in 2:nrow(C)) {
      new.mat_group[, , c_count] <- t((new.mat[c((e + 1):(c_count * ncol(y))), ]))
      beta.group[, , c_count] <- t(Big_beta_respone[c((e + 1):(c_count * ncol(y))), ])

      beta_transform <- matrix(0, p, (K + 1) * ncol(y))
      beta_transform[, 1:(1 + K)] <- matrix(new.mat_group[, 1, c_count], ncol = (K + 1), nrow = p)
      input2 <- 1:(ncol(y) * (1 + K))
      multiple_of_K <- (input2 %% (K + 1)) == 0
      II2 <- input2[multiple_of_K]
      e2 <- II2[-length(II2)][1]

      for (c_count2 in 2:ncol(y)) {
        beta_transform[, c((e2 + 1):(c_count2 * (1 + K)))] <- matrix(new.mat_group[, c_count2, c_count], ncol = (K + 1), nrow = p)
        e2 <- II2[c_count2]
      }

      norm_res <- ((apply(beta_transform, 1, twonorm)))
      coef.term1 <- pmax(1 - (gg[1]) / rho / (norm_res), 0)

      N_E1 <- scale(t(beta_transform), center = FALSE, scale = 1 / coef.term1)

      N_E1 <- t(N_E1)
      beta_transform1 <- matrix(0, p + p * K, ncol(y))
      beta_transform1[, 1] <- as.vector(N_E1[, 1:(K + 1)])

      input3 <- 1:(ncol(y) * (1 + K))
      multiple_of_K <- (input3 %% (K + 1)) == 0
      II3 <- input3[multiple_of_K]
      e3 <- II3[-length(II3)][1]

      for (c_count3 in 2:ncol(y)) {
        beta_transform1[, c_count3] <- as.vector(N_E1[, c((e3 + 1):((K + 1) * c_count3))])
        e3 <- II3[c_count3]
      }

      N_E[[c_count]] <- (t((beta_transform1)))

      e <- II[c_count]
    }

    E[seq_len(ncol(C)), ] <- N_E[[1]]

    c_count <- 2
    e <- II[-length(II)][1]
    for (c_count in 2:nrow(C)) {
      E[c((e + 1):(c_count * ncol(y))), ] <- N_E[[c_count]]
      e <- II[c_count]
    }

    H <- H + Big_beta_respone - E

    obj <- c(obj, obj)

    # Calculate residuals for iteration t
    v.diff <- sum((-rho * (V - V_old))^2, na.rm = TRUE)

    q.diff <- sum((-rho * (Q - Q_old))^2, na.rm = TRUE)
    e.diff <- sum((-rho * (E - E_old))^2, na.rm = TRUE)
    ee.diff <- sum((-rho * (EE - EE_old))^2, na.rm = TRUE)

    s <- sqrt(v.diff + q.diff + e.diff + ee.diff)

    v.diff1 <- sum(v.diff1)
    q.diff1 <- sum(q.diff1)
    e.diff1 <- sum(((Big_beta_respone - E))^2, na.rm = TRUE)
    ee.diff1 <- sum(ee.diff1)
    r <- sqrt(v.diff1 + q.diff1 + e.diff1 + ee.diff1)

    res_dual <- s
    res_pri <- r

    e.primal <- sqrt(length(Big_beta11) + 2 * length(beta_hat) + length(Big_beta_respone)) * e.abs + e.rel * max(twonorm(c((Big_beta11), (beta_hat), (beta_hat), (Big_beta_respone))), twonorm(-c((V), (Q), (E), (EE))))

    e.dual <- sqrt(length(Big_beta11) + 2 * length(beta_hat) + length(Big_beta_respone)) * e.abs + e.rel * twonorm((c((O), (P), (H), (HH))))
    V_old <- V
    Q_old <- Q
    E_old <- E
    EE_old <- EE

    if (res_pri > 10 * res_dual) {
      rho <- 2 * rho
    } else if (res_pri * 10 < res_dual) {
      rho <- rho / 2
    }

    if (my_print) {
      cat(res_dual, e.dual, res_pri, e.primal)
    }
    if (res_pri <= e.primal && res_dual <= e.dual) {
      # Remove excess beta and nll

      # Update convergence message
      message("Convergence reached after ", i, " iterations")
      converge <- TRUE
      break
    }
    converge <- FALSE
  } ### iteration

  res_val <- t(I) %*% (E)
  for (jj in seq_len(ncol(y))) {
    group <- (t(G) %*% t((V[, , jj])))

    group1 <- group[1, ]
    group2 <- t(group[-1, ])
    new_group <- matrix(0, p, (K + 1))
    new_group[, 1] <- group1
    new_group[, -1] <- group2
    new_g_theta <- as.vector(new_group)

    finB1 <- as.vector(beta_hat[1:p, jj]) * (new_g_theta[1:p] != 0) * (as.vector((Q[, 1, jj])) != 0)
    finB2 <- as.vector(beta_hat[-1:-p, jj]) * (new_g_theta[-1:-p] != 0) * (as.vector((Q[, -1, jj])) != 0)

    beta_hat1 <- matrix(c(finB1, finB2), ncol = (K + 1), nrow = p)
    beta[, jj] <- beta_hat1[, 1]
    theta[, , jj] <- (beta_hat1[, -1])
    beta_hat[, jj] <- c(c(beta_hat1[, 1], as.vector(theta[, , jj])))
  }
  y_hat <- model_p(beta0, theta0, beta = beta_hat, X = W_hat, Z)

  out <- list(beta0 = beta0, theta0 = theta0, beta = beta, theta = theta, converge = converge, obj = obj, beta_hat = beta_hat, y_hat = y_hat)

  return(out)
}
