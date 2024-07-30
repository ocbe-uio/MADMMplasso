hh_nlambda_loop <- function(
  lam, nlambda, beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it,
  my_W_hat, XtY, y, N, e.abs, e.rel, alpha, alph, svd.w, tree, my_print,
  invmat, gg, tol, parallel, pal, BETA0, THETA0, BETA,
  BETA_hat, Y_HAT, THETA, D, my_values
) {
  obj <- NULL
  non_zero_theta <- NULL
  my_obj <- list()
  n_main_terms <- NULL
  lam_list <- list()
  hh <- 1
  while (hh <= nlambda) {
    lambda <- lam[hh, ]

    start_time <- Sys.time()
    if (pal) {
      my_values <- admm_MADMMplasso(
        beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
        y, N, e.abs, e.rel, alpha, lambda, alph, svd.w, tree, my_print, invmat,
        gg[hh, ]
      )
      beta <- my_values$beta
      theta <- my_values$theta
      my_obj[[hh]] <- list(my_values$obj)
      beta0 <- my_values$beta0
      theta0 <- my_values$theta0 ### iteration
      beta_hat <- my_values$beta_hat
      y_hat <- my_values$y_hat
    }
    if (my_print) {
      print(Sys.time() - start_time)
    }
    if (parallel && !pal) {
      beta <- my_values[[hh]]$beta
      theta <- my_values[[hh]]$theta
      my_obj[[hh]] <- list(my_values[[hh]]$obj)
      beta0 <- my_values[[hh]]$beta0
      theta0 <- my_values[[hh]]$theta0 ### iteration
      beta_hat <- my_values[[hh]]$beta_hat
      y_hat <- my_values[[hh]]$y_hat
    } else if (!parallel && !pal) {
      beta <- my_values[[hh]]$beta
      theta <- my_values[[hh]]$theta
      my_obj[[hh]] <- list(my_values[[hh]]$obj)
      beta0 <- my_values[[hh]]$beta0
      theta0 <- my_values[[hh]]$theta0 ### iteration
      beta_hat <- my_values[[hh]]$beta_hat
      y_hat <- my_values[[hh]]$y_hat
    }
    # Executed if par == TRUE, independent of parallel

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

    if (my_print) {
      if (hh == 1) {
        cat(hh, (n_main_terms[hh]), non_zero_theta[hh], obj1)
      } else {
        cat(hh, (n_main_terms[hh]), non_zero_theta[hh], obj[hh - 1], obj1)
      }
    }

    hh <- hh + 1
  } ### lambda

  # Output
  list(
    obj = obj, n_main_terms = n_main_terms, non_zero_theta = non_zero_theta,
    BETA0 = BETA0, THETA0 = THETA0, BETA = BETA, BETA_hat = BETA_hat,
    Y_HAT = Y_HAT, THETA = THETA
  )
}
