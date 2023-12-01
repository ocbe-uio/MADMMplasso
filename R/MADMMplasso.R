#' @title Fit a multi-response pliable lasso model over a path of regularization values
#' @description TODO: add description (This function fits a multi-response pliable lasso model over a path of regularization values?)
#' @param X  N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical varables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param y N by D matrix  of responses. The X and Z variables are centered in the function. We recommmend that X and Z also be standardized before the call
#' @param nlambda number of lambda values desired (default 50).
#' @param alpha mixing parameter- default 0.5
#' @param max_it maximum number of iterations in the ADMM algorithm for one lambda. Default 50000
#' @param maxgrid similar to nlambda
#' @param rho the Lagrange variable for the ADMM
#' @param e.abs absolute error for the admm
#' @param e.rel relative error for the admm
#' @param gg penalty term for the tree structure
#' @param my_lambda TODO: fill in
#' @param lambda_min TODO: fill in
#' @param max_it TODO: fill in
#' @param my_print TODO: fill in
#' @param alph TODO: fill in
#' @param tree TODO: fill in
#' @param cv TODO: fill in
#' @param parallel TODO: fill in
#' @param pal TODO: fill in
#' @param tol TODO: fill in
#' @param cl TODO: fill in
#' @return  predicted values for the MADMMplasso fit
#' @example inst/examples/MADMMplasso_example.R
#' @export
MADMMplasso <- function(X, Z, y, alpha, my_lambda = NULL, lambda_min = .001, max_it = 50000, e.abs = 1E-3, e.rel = 1E-3, maxgrid, nlambda, rho = 5, my_print = F, alph = 1.8, tree, cv = F, parallel = T, pal = 0, gg = NULL, tol = 1E-4, cl = 4,legacy=F) {
  N <- nrow(X)

  p <- ncol(X)
  K <- ncol(Z)
  D <- dim(y)[2]

  TT <- tree

  C <- TT$Tree
  CW <- TT$Tw

  BETA0 <- lapply(
    seq_len(nlambda),
    function(j) (matrix(0, nrow = (D)))
  )

  BETA <- lapply(
    seq_len(nlambda),
    function(j) (as(matrix(0, nrow = p, ncol = (D)), "sparseMatrix"))
  )
  BETA_hat <- lapply(
    seq_len(nlambda),
    function(j) (as(matrix(0, nrow = p + p * K, ncol = (D)), "sparseMatrix"))
  )

  THETA0 <- lapply(
    seq_len(nlambda),
    function(j) (matrix(0, nrow = K, ncol = (D)))
  )

  THETA <- lapply(
    seq_len(nlambda),
    function(j) (as.sparse3Darray(array(0, c(p, K, (D)))))
  )

  Y_HAT <- lapply(
    seq_len(nlambda),
    function(j) (matrix(0, nrow = N, ncol = (D)))
  )
  rat <- lambda_min

  if (is.null(my_lambda)) {
    lamda_new <- matrix(0, dim(y)[2])
    r <- y

    lammax <- lapply(
      seq_len(dim(y)[2]),
      function(g) {
        l_max <- max(abs(t(X) %*% (r - colMeans(r))) / length(r[, 1])) / ((1 - alpha) + (max(gg[1, ]) * max(CW) + max(gg[2, ])))
        return(l_max)
      }
    )

    big_lambda <- lammax
    lambda_i <- lapply(
      seq_len(dim(y)[2]),
      function(g) {
        lam_i <- exp(seq(log(big_lambda[[g]]), log(big_lambda[[g]] * rat), length = maxgrid))

        return(lam_i)
      }
    )
    gg1 <- seq((gg[1, 1]), (gg[1, 2]), length = maxgrid)
    gg2 <- seq((gg[2, 1]), (gg[2, 2]), length = maxgrid)
    gg3 <- matrix(0, maxgrid, 2)
    gg3[, 1] <- gg1
    gg3[, 2] <- gg2
    gg <- gg3
  } else {
    lambda_i <- my_lambda
    gg1 <- gg
    gg <- gg1
  }

  lam_list <- list()
  beta_0_list <- list()
  theta_0_list <- list()
  beta_list <- list()
  theta_list <- list()
  obj <- c()
  n_main_terms <- c()
  non_zero_theta <- c()
  my_obj <- list()

  my_W_hat <- generate_my_w(X = X, Z = Z, quad = TRUE)

  svd.w <- svd(my_W_hat)
  svd.w$tu <- t(svd.w$u)
  svd.w$tv <- t(svd.w$v)

  rho1 <- rho

  D <- dim(y)[2]

  ### for response groups ###############################################################
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
    DD1 <- rho1 * (new_I[rr] + 1)
    DD2 <- new_G + DD1
    invmat[[rr]] <- DD2
  }

  beta0 <- matrix(0, 1, D)
  theta0 <- matrix(0, K, D)
  beta <- (matrix(0, p, D))
  beta_hat <- (matrix(0, p + p * (K), D))

  # auxiliary variables for the L1 norm####

  theta <- (array(0, c(p, K, D)))

  if (is.null(my_lambda)) {
    lam <- matrix(0, nlambda, dim(y)[2])
    for (i in 1:dim(y)[2]) {
      lam[, i] <- lambda_i[[i]]
    }
  } else {
    lam <- lambda_i
  }

  r_current <- y
  b <- reg(r_current, Z)
  beta0 <- b$beta0
  theta0 <- b$theta0

  new_y <- y - (matrix(1, N) %*% beta0 + Z %*% ((theta0)))

  XtY <- crossprod((my_W_hat), (new_y))

  cl1 <- cl
  if (parallel) {
    cl <- makeCluster(cl1, type = "FORK")

    doParallel::registerDoParallel(cl = cl)
    foreach::getDoParRegistered()

    my_values <- foreach(i = 1:nlambda, .packages = "MADMMplasso", .combine = rbind) %dopar% {
      admm_MADMMplasso(
        beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
        y, N, e.abs, e.rel, alpha, lam[i, ], alph, svd.w, tree, my_print,
        invmat, cv, gg[i, ],legacy=legacy
      )
    }
    parallel::stopCluster(cl)
  } else if (parallel == F & pal == 0) {
    my_values <- lapply(
      seq_len(nlambda),
      function(g) {
        admm_MADMMplasso(
          beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat,
          XtY, y, N, e.abs, e.rel, alpha, lam[g, ], alph, svd.w, tree, my_print,
          invmat, cv, gg[g, ],legacy=legacy
        )
      }
    )
  }

  repeat_loop <- 0
  hh <- 1
  while (hh <= nlambda) {
    res_dual <- 0 # dual residual
    res_pri <- 0 # primal residual

    lambda <- lam[hh, ]

    start_time <- Sys.time()
    if (pal == 1) {
      my_values <- admm_MADMMplasso(
        beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
        y, N, e.abs, e.rel, alpha, lambda, alph, svd.w, tree, my_print, invmat,
        cv, gg[hh, ],legacy=legacy
      )

      beta <- my_values$beta
      theta <- my_values$theta
      converge <- my_values$converge
      my_obj[[hh]] <- list(my_values$obj)
      beta0 <- my_values$beta0
      theta0 <- my_values$theta0 ### iteration
      beta_hat <- my_values$beta_hat
      y_hat <- my_values$y_hat
    }
    cost_time <- Sys.time() - start_time
    print(cost_time)
    if (parallel == T & pal == 0) {
      beta <- my_values[hh, ]$beta
      theta <- my_values[hh, ]$theta
      converge <- my_values[hh, ]$converge
      my_obj[[hh]] <- list(my_values[hh, ]$obj)
      beta0 <- my_values[hh, ]$beta0
      theta0 <- my_values[hh, ]$theta0 ### iteration
      beta_hat <- my_values[hh, ]$beta_hat
      y_hat <- my_values[hh, ]$y_hat
    } else if (parallel == F & pal == 0) {
      beta <- my_values[[hh]]$beta
      theta <- my_values[[hh]]$theta
      converge <- my_values[[hh]]$converge
      my_obj[[hh]] <- list(my_values[[hh]]$obj)
      beta0 <- my_values[[hh]]$beta0
      theta0 <- my_values[[hh]]$theta0 ### iteration
      beta_hat <- my_values[[hh]]$beta_hat
      y_hat <- my_values[[hh]]$y_hat
    }

    beta1 <- as(beta * (abs(beta) > tol), "sparseMatrix")
    theta1 <- as.sparse3Darray(theta * (abs(theta) > tol))
    beta_hat1 <- as(beta_hat * (abs(beta_hat) > tol), "sparseMatrix")

    n_interaction_terms <- count_nonzero_a((theta1))

    n_main_terms <- (c(n_main_terms, count_nonzero_a((beta1))))

    obj1 <- (sum(as.vector((y - y_hat)^2))) / (D * N)
    obj <- c(obj, obj1)

    non_zero_theta <- (c(non_zero_theta, n_interaction_terms))
    lam_list <- (c(lam_list, lambda))

    BETA0[[hh]] <- beta0
    THETA0[[hh]] <- theta0
    BETA[[hh]] <- as(beta1, "sparseMatrix")
    BETA_hat[[hh]] <- as(beta_hat1, "sparseMatrix")

    Y_HAT[[hh]] <- y_hat
    THETA[[hh]] <- as.sparse3Darray(theta1)

    if (cv == F) {
      if (hh == 1) {
        print(c(hh, (n_main_terms[hh]), non_zero_theta[hh], obj1))
      } else {
        print(c(hh, (n_main_terms[hh]), non_zero_theta[hh], obj[hh - 1], obj1))
      }
    } else {
      if (hh == 1) {
        print(c(hh, (n_main_terms[hh]), non_zero_theta[hh], obj1))
      } else {
        print(c(hh, (n_main_terms[hh]), non_zero_theta[hh], obj[hh - 1], obj1))
      }
    }

    hh <- hh + 1
  } ### lambda

  remove(invmat)
  remove(my_values)
  remove(my_W_hat)

  obj[1] <- obj[2]

  pred <- data.frame(Lambda = lam, nzero = n_main_terms, nzero_inter = non_zero_theta, OBJ_main = obj)
  out <- list(beta0 = BETA0, beta = BETA, BETA_hat = BETA_hat, theta0 = THETA0, theta = THETA, path = pred, Lambdas = lam, non_zero = n_main_terms, LOSS = obj, it.obj = my_obj, Y_HAT = Y_HAT, gg = gg)
  class(out) <- "MADMMplasso"
  # Return results
  return(out)
}
