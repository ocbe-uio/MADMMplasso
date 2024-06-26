#include "MADMMplasso.h"
#include "misc.h"

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
  const bool my_print,
  const arma::mat gg,
  const double tol,
  const bool parallel,
  const bool pal,
  arma::cube BETA0,
  arma::cube THETA0,
  arma::cube BETA,
  arma::cube BETA_hat,
  arma::cube Y_HAT,
  const unsigned int D,
  const arma::sp_mat C,
  const arma::vec CW,
  const arma::mat svd_w_tu,
  const arma::mat svd_w_tv,
  const arma::vec svd_w_d,
  Rcpp::List my_values
) {
  arma::vec obj;
  arma::vec non_zero_theta;
  arma::vec n_main_terms;
  arma::vec lam_list;
  arma::mat y_hat = y;
  unsigned int hh = 0;
  arma::field<arma::cube> THETA(nlambda);
  while (hh <= nlambda - 1) {
    arma::vec lambda = lam.row(hh).t();

    if (pal) {
      // In this case, my_values is an empty list to be created now
      arma::field<arma::cube> my_values_hh = admm_MADMMplasso_cpp(
        beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, my_W_hat, XtY,
        y, N, e_abs, e_rel, alpha, lambda, alph, svd_w_tu, svd_w_tv, svd_w_d, C, CW,
        gg.row(hh), my_print
      );
      beta0 = my_values_hh(0).slice(0);
      theta0 = my_values_hh(1).slice(0);
      beta = my_values_hh(2).slice(0);
      theta = my_values_hh(3);
      beta_hat = my_values_hh(5).slice(0);
      y_hat = my_values_hh(6).slice(0);
    } else {
      // Gets triggered regardless of parallel. Whatever the case,
      // my_values is already a list of length hh
      arma::field<arma::cube> my_values_hh = my_values[hh];
      beta0 = my_values_hh(0).slice(0);
      theta0 = my_values_hh(1).slice(0);
      beta = my_values_hh(2).slice(0);
      theta = my_values_hh(3);
      beta_hat = my_values_hh(5).slice(0);
      y_hat = my_values_hh(6).slice(0);
    }

    // should be sparse, but Arma doesn't have sp_cube; beta1 and beta_hat1
    // are going into a cube, so they need to be dense as well
    arma::mat beta1(beta % (abs(beta) > tol));
    arma::cube theta1(theta % (abs(theta) > tol));
    arma::mat beta_hat1(beta_hat % (abs(beta_hat) > tol));

    // TODO: messy! Simplify!
    arma::vec n_interaction_terms(1);
    n_interaction_terms = count_nonzero_a_cube(theta1);
    arma::vec n_beta_terms(1);
    n_beta_terms = count_nonzero_a_mat(beta1);
    n_main_terms = arma::join_vert(n_main_terms, n_beta_terms);

    double obj1 = arma::accu(arma::pow(y - y_hat, 2)) / (D * N);
    obj.resize(obj.n_elem + 1);
    obj(obj.n_elem - 1) = obj1;
    non_zero_theta = arma::join_vert(non_zero_theta, n_interaction_terms);
    lam_list = arma::join_vert(lam_list, lambda);

    BETA0.slice(hh) = beta0;
    THETA0.slice(hh) = theta0;
    BETA.slice(hh) = beta1;
    BETA_hat.slice(hh) = beta_hat1;
    Y_HAT.slice(hh) = y_hat;
    THETA(hh) = theta1;

    if (my_print) {
      if (hh == 0) {
        Rcpp::Rcout << hh << n_main_terms(hh) << non_zero_theta(hh) << obj1 << std::endl;
      } else {
        Rcpp::Rcout << hh << n_main_terms(hh) << non_zero_theta(hh) << obj(hh - 1) << obj1 << std::endl;
      }
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
