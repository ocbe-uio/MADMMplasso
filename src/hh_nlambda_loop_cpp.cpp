#include "MADMMplasso.h"

// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::List hh_nlambda_loop_cpp(
  arma::mat lam,
  unsigned int nlambda,
  const arma::vec beta0,
  const arma::mat theta0,
  arma::mat beta,
  arma::mat beta_hat,
  arma::cube theta,
  const double rho1,
  const arma::mat X,
  const arma::mat Z,
  const int max_it,
  const arma::mat my_W_hat,
  const arma::mat XtY,
  const arma::mat y,
  const int N,
  const double e_abs,
  const double e_rel,
  const double alpha,
  const double alph,
  const Rcpp::List svd_w,
  const Rcpp::List tree,
  const bool my_print,
  const Rcpp::List invmat,
  const arma::mat gg,
  const double tol,
  const bool parallel,
  const bool pal,
  Rcpp::List BETA0,
  Rcpp::List THETA0,
  Rcpp::List BETA,
  Rcpp::List BETA_hat,
  Rcpp::List Y_HAT,
  Rcpp::List THETA,
  const unsigned int D,
  Rcpp::List my_values
) {
  arma::vec obj;
  arma::vec non_zero_theta;
  Rcpp::List my_obj;
  arma::vec n_main_terms;
  Rcpp::List lam_list;
  unsigned int hh = 1;
  while (hh <= nlambda) {
    arma::vec lambda = lam.row(hh);

    // start_time <- Sys.time()
    if (pal) {
      // my_values <- admm_MADMMplasso(
      //   beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
      //   y, N, e_abs, e_rel, alpha, lambda, alph, svd_w, tree, my_print, invmat,
      //   gg[hh, ], legacy
      // )

      // beta <- my_values$beta
      // theta <- my_values$theta
      // my_obj[[hh]] <- list(my_values$obj)
      // beta0 <- my_values$beta0
      // theta0 <- my_values$theta0 ### iteration
      // beta_hat <- my_values$beta_hat
      // y_hat <- my_values$y_hat
    } else if (parallel) {
      // beta <- my_values[hh, ]$beta
      // theta <- my_values[hh, ]$theta
      // my_obj[[hh]] <- list(my_values[hh, ]$obj)
      // beta0 <- my_values[hh, ]$beta0
      // theta0 <- my_values[hh, ]$theta0 ### iteration
      // beta_hat <- my_values[hh, ]$beta_hat
      // y_hat <- my_values[hh, ]$y_hat
    } else {
      // beta <- my_values[[hh]]$beta
      // theta <- my_values[[hh]]$theta
      // my_obj[[hh]] <- list(my_values[[hh]]$obj)
      // beta0 <- my_values[[hh]]$beta0
      // theta0 <- my_values[[hh]]$theta0 ### iteration
      // beta_hat <- my_values[[hh]]$beta_hat
      // y_hat <- my_values[[hh]]$y_hat
    }

    // beta1 <- as(beta * (abs(beta) > tol), "sparseMatrix")
    // theta1 <- as.sparse3Darray(theta * (abs(theta) > tol))
    // beta_hat1 <- as(beta_hat * (abs(beta_hat) > tol), "sparseMatrix")

    // n_interaction_terms <- count_nonzero_a((theta1))

    // n_main_terms <- (c(n_main_terms, count_nonzero_a((beta1))))

    // obj1 <- (sum(as.vector((y - y_hat)^2))) / (D * N)
    double obj1 = 0; // TEMP
    // obj <- c(obj, obj1)

    // non_zero_theta <- (c(non_zero_theta, n_interaction_terms))
    // lam_list <- (c(lam_list, lambda))

    // BETA0[[hh]] <- beta0
    // THETA0[[hh]] <- theta0
    // BETA[[hh]] <- as(beta1, "sparseMatrix")
    // BETA_hat[[hh]] <- as(beta_hat1, "sparseMatrix")

    // Y_HAT[[hh]] <- y_hat
    // THETA[[hh]] <- as.sparse3Darray(theta1)

    if (hh == 1) {
      Rcpp::Rcout << hh << n_main_terms[hh] << non_zero_theta[hh] << obj1 << std::endl;
    } else {
      Rcpp::Rcout << hh << n_main_terms[hh] << non_zero_theta[hh] << obj[hh - 1] << obj1 << std::endl;
    }
    hh += 1;
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("obj") = obj,
    Rcpp::Named("n_main_terms") = n_main_terms,
    Rcpp::Named("non_zero_theta") = non_zero_theta,
    Rcpp::Named("BETA0") = BETA0,
    Rcpp::Named("THETA0") = THETA0,
    Rcpp::Named("BETA") = BETA,
    Rcpp::Named("BETA_hat") = BETA_hat,
    Rcpp::Named("Y_HAT") = Y_HAT,
    Rcpp::Named("THETA") = THETA
  );
  return out;
}
