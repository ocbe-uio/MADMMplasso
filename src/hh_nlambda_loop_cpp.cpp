#include "MADMMplasso.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List hh_nlambda_loop_cpp(
  const arma::mat lam,
  const unsigned int nlambda,
  arma::vec beta0,
  arma::mat theta0,
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
  arma::vec lam_list;
  arma::mat y_hat = y;
  unsigned int hh = 0;

  Rcpp::Environment MADMMplasso = Rcpp::Environment::namespace_env("MADMMplasso"); // TEMP
  Rcpp::Function count_nonzero_a = MADMMplasso["count_nonzero_a"]; // TEMP
  while (hh <= nlambda - 1) {
    arma::vec lambda = lam.row(hh).t();

    Rcpp::List my_values_hh;
    if (parallel) {
      // my_values is already a list of length hh
      my_values_hh = my_values[hh];
    } else if (pal) {
      // In this case, my_values is an empty list to be created now
      my_values_hh = admm_MADMMplasso_cpp(
        beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
        y, N, e_abs, e_rel, alpha, lambda, alph, svd_w, tree, invmat,
        gg.row(hh).t(), my_print
      );
    }

    beta = Rcpp::as<arma::mat>(my_values_hh["beta"]);
    theta = Rcpp::as<arma::cube>(my_values_hh["theta"]);
    beta0 = Rcpp::as<arma::vec>(my_values_hh["beta0"]);
    theta0 = Rcpp::as<arma::mat>(my_values_hh["theta0"]);
    beta_hat = Rcpp::as<arma::mat>(my_values_hh["beta_hat"]);
    y_hat = Rcpp::as<arma::mat>(my_values_hh["y_hat"]);

    arma::sp_mat beta1(beta % (abs(beta) > tol));
    arma::cube theta1(theta % (abs(theta) > tol)); // should be sparse, but Arma doesn't have sp_cube
    arma::sp_mat beta_hat1(beta_hat % (abs(beta_hat) > tol));

    arma::vec n_interaction_terms = Rcpp::as<arma::vec>(count_nonzero_a(theta1));
    n_main_terms = arma::join_vert(n_main_terms, Rcpp::as<arma::vec>(count_nonzero_a(beta1)));

    double obj1 = arma::accu(arma::pow(y - y_hat, 2)) / (D * N);
    obj.resize(obj.n_elem + 1);
    obj(obj.n_elem - 1) = obj1;
    non_zero_theta = arma::join_vert(non_zero_theta, n_interaction_terms);
    lam_list = arma::join_vert(lam_list, lambda);

    BETA0[hh] = beta0;
    THETA0[hh] = theta0;
    BETA[hh] = arma::conv_to<arma::sp_mat>::from(beta1);
    BETA_hat[hh] = arma::conv_to<arma::sp_mat>::from(beta_hat1);
    Y_HAT[hh] = y_hat;
    THETA[hh] = theta1;

    if (hh == 0) {
      Rcpp::Rcout << hh << n_main_terms(hh) << non_zero_theta(hh) << obj1 << std::endl;
    } else {
      Rcpp::Rcout << hh << n_main_terms(hh) << non_zero_theta(hh) << obj(hh - 1) << obj1 << std::endl;
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
