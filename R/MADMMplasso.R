
#' @title Fit a multi-response pliable lasso model over a path of regularization values
#' @description This function fits a multi-response pliable lasso model over a path of regularization values.
#' @param X  N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical variables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param y N by D matrix  of responses. The X and Z variables are centered in the function. We recommend that X and Z also be standardized before the call
#' @param maxgrid  number of lambda_3 values desired
#' @param nlambda  number of lambda_3 values desired. Similar to maxgrid but can have a value less than or equal to maxgrid.
#' @param alpha mixing parameter. When the goal is to include more interactions, alpha should be very small and vice versa.
#' @param max_it maximum number of iterations in the ADMM algorithm for one lambda
#' @param rho the Lagrange variable for the ADMM. This value is updated during the ADMM call based on a certain condition.
#' @param e.abs absolute error for the ADMM
#' @param e.rel relative error for the ADMM
#' @param gg penalty term for the tree structure. This is a 2×2 matrix values in the first row representing the maximum to the minimum values for lambda_1 and the second row representing the maximum to the minimum values for lambda_2. In the current setting, we set both maximum and the minimum to be same because cross validation is not carried across the lambda_1 and lambda_2. However, setting different values will work during the model fit.
#' @param my_lambda user specified lambda_3 values
#' @param lambda_min the smallest value for lambda_3 , as a fraction of max(lambda_3), the (data derived (lammax)) entry value (i.e. the smallest value for which all coefficients are zero)
#' @param max_it 	maximum number of iterations in loop for one lambda during the ADMM optimization
#' @param my_print Should information form each ADMM iteration be printed along the way? This prints the dual and primal residuals
#' @param alph an overrelaxation parameter in \[1, 1.8\]. The implementation is borrowed from Stephen Boyd's \href{https://stanford.edu/~boyd/papers/admm/lasso/lasso.html}{MATLAB code}
#' @param tree The results from the hierarchical clustering of the response matrix. The easy way to obtain this is by using the function (tree_parms) which gives a default clustering. However, user decide on a specific structure and then input a tree that follows such structure.
#' @param pal Should the lapply function be applied for an alternative to parallelization.
#' @param tol threshold for the non-zero coefficients
#' @param cl The number of CPUs to be used for parallel processing
#' @param legacy If \code{TRUE}, use the R version of the algorithm
#' @return   predicted values for the MADMMplasso object with the following components:
#' path: a table containing the summary of the model for each lambda_3.
#'
#' beta0: a list (length=nlambda) of estimated beta_0 coefficients each having a size of 1 by ncol(y)
#'
#' beta: a list (length=nlambda) of estimated beta coefficients each having a matrix   ncol(X) by ncol(y)
#'
#' BETA_hat: a list (length=nlambda) of estimated beta and theta coefficients each having a matrix   (ncol(X)+ncol(X) by ncol(Z)) by ncol(y)
#'
#'  theta0: a list (length=nlambda) of estimated theta_0 coefficients each having a matrix   ncol(Z) by ncol(y)
#'
#'  theta: a list (length=nlambda) of estimated theta coefficients each having a an array   ncol(X) by ncol(Z) by ncol(y)
#'
#'  Lambdas: values of lambda_3 used
#'
#'  non_zero: number of nonzero betas for each model in path
#'
#'  LOSS: sum of squared of the error for each model in path
#'
#'   Y_HAT: a list (length=nlambda) of predicted response nrow(X) by ncol(y)
#'
#'  gg: penalty term for the tree structure (lambda_1 and lambda_2) for each lambda_3 nlambda by 2

#' @example inst/examples/MADMMplasso_example.R
#' @export
MADMMplasso <- function(X, Z, y, alpha, my_lambda = NULL, lambda_min = 0.001, max_it = 50000, e.abs = 1E-3, e.rel = 1E-3, maxgrid, nlambda, rho = 5, my_print = FALSE, alph = 1.8, tree, pal = cl == 1L, gg = NULL, tol = 1E-4, cl = 1L, legacy = FALSE) {
  # Recalculating the number of CPUs
  if (pal && cl > 1L) {
    cl <- 1L
    warning("pal is TRUE, resetting cl to 1")
  }
  parallel <- cl > 1L
  if (my_print) {
    message(
      "Parallelization is ", ifelse(parallel, "enabled ", "disabled "),
      "(", cl, " CPUs)"
    )
    message("pal is ", ifelse(pal, "TRUE", "FALSE"))
    message("Using ", ifelse(legacy, "R", "C++"), " engine")
  }
  N <- nrow(X)

  p <- ncol(X)
  K <- ncol(Z)
  D <- ncol(y)

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
    # Validation
    stopifnot(nlambda <= maxgrid)

    r <- y

    lammax <- lapply(
      seq_len(ncol(y)),
      function(g) {
        max(abs(t(X) %*% (r - colMeans(r))) / length(r[, 1])) / ((1 - alpha) + (max(gg[1, ]) * max(CW) + max(gg[2, ])))
      }
    )

    big_lambda <- lammax
    lambda_i <- lapply(
      seq_len(ncol(y)),
      function(g) {
        exp(seq(log(big_lambda[[g]]), log(big_lambda[[g]] * rat), length = maxgrid))
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

  my_W_hat <- generate_my_w(X = X, Z = Z)

  svd.w <- svd(my_W_hat)
  svd.w$tu <- t(svd.w$u)
  svd.w$tv <- t(svd.w$v)

  rho1 <- rho

  D <- ncol(y)

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
  new_G <- matrix(0, (p + p * K))
  new_G[1:p] <- 1
  new_G[-1:-p] <- 2
  new_G[1:p] <- rho * (1 + new_G[1:p])
  new_G[-1:-p] <- rho * (1 + new_G[-1:-p])

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
    lam <- matrix(0, nlambda, ncol(y))
    for (i in seq_len(ncol(y))) {
      lam[, i] <- lambda_i[[i]]
    }
  } else {
    lam <- lambda_i
  }

  r_current <- y
  b <- reg(r_current, Z)
  beta0 <- b[1, ]
  theta0 <- b[-1, ]

  new_y <- y - (matrix(1, N) %*% beta0 + Z %*% ((theta0)))

  XtY <- crossprod((my_W_hat), (new_y))


  cl1 <- cl

  # Adjusting objects for C++
  if (!legacy) {
    C <- TT$Tree
    CW <- TT$Tw
    svd_w_tu <- t(svd.w$u)
    svd_w_tv <- t(svd.w$v)
    svd_w_d <- svd.w$d
    BETA <- array(0, c(p, D, nlambda))
    BETA_hat <- array(0, c(p + p * K, D, nlambda))
  }

  # Pre-calculating my_values through my_values_matrix
  if (parallel) {
    if (.Platform$OS.type == "unix") {
      cl <- parallel::makeForkCluster(cl1)
    } else {
      cl <- parallel::makeCluster(cl1)
    }
    doParallel::registerDoParallel(cl = cl)
    foreach::getDoParRegistered()
    if (legacy) {
      my_values <- foreach(i = 1:nlambda) %dopar% {
        admm_MADMMplasso(
          beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
          y, N, e.abs, e.rel, alpha, lam[i, ], alph, svd.w, tree, my_print,
          invmat, gg[i, ]
        )
      }
    } else {
      my_values <- foreach(i = 1:nlambda) %dopar% {
        admm_MADMMplasso_cpp(
          beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
          y, N, e.abs, e.rel, alpha, lam[i, ], alph, svd_w_tu, svd_w_tv, svd_w_d,
          C, CW, gg[i, ], my_print
        )
      }
    }
    parallel::stopCluster(cl)
  } else if (!parallel && !pal) {
    if (legacy) {
      my_values <- lapply(
        seq_len(nlambda),
        function(g) {
          admm_MADMMplasso(
            beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat,
            XtY, y, N, e.abs, e.rel, alpha, lam[g, ], alph, svd.w, tree, my_print,
            invmat, gg[g, ]
          )
        }
      )
    } else {
      my_values <- lapply(
        seq_len(nlambda),
        function(i) {
          admm_MADMMplasso_cpp(
            beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
            y, N, e.abs, e.rel, alpha, lam[i, ], alph, svd_w_tu, svd_w_tv, svd_w_d,
            C, CW, gg[i, ], my_print
          )
        }
      )
    }
  } else {
    # This is triggered when parallel is FALSE and pal is TRUE
    my_values <- list()
  }

  # Big calculations
  if (legacy) {
    loop_output <- hh_nlambda_loop(
      lam, nlambda, beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it,
      my_W_hat, XtY, y, N, e.abs, e.rel, alpha, alph, svd.w, tree, my_print,
      invmat, gg, tol, parallel, pal, BETA0, THETA0, BETA,
      BETA_hat, Y_HAT, THETA, D, my_values
    )
  } else {
    loop_output <- hh_nlambda_loop_cpp(
      lam, as.integer(nlambda), beta0, theta0, beta, beta_hat, theta, rho1, X, Z, as.integer(max_it),
      my_W_hat, XtY, y, as.integer(N), e.abs, e.rel, alpha, alph, my_print,
      gg, tol, parallel, pal, simplify2array(BETA0), simplify2array(THETA0),
      BETA, BETA_hat, simplify2array(Y_HAT),
      as.integer(D), C, CW, svd_w_tu, svd_w_tv, svd_w_d, my_values
    )
    loop_output <- post_process_cpp(loop_output)
  }

  # Final adjustments in output
  loop_output$obj[1] <- loop_output$obj[2]

  pred <- data.frame(
    Lambda = lam,
    nzero = loop_output$n_main_terms,
    nzero_inter = loop_output$non_zero_theta,
    OBJ_main = loop_output$obj
  )
  out <- list(
    beta0 = loop_output$BETA0,
    beta = loop_output$BETA,
    BETA_hat = loop_output$BETA_hat,
    theta0 = loop_output$THETA0,
    theta = loop_output$THETA,
    path = pred,
    Lambdas = lam,
    non_zero = loop_output$n_main_terms,
    LOSS = loop_output$obj,
    Y_HAT = loop_output$Y_HAT,
    gg = gg
  )
  class(out) <- "MADMMplasso"
  # Return results
  return(out)
}
