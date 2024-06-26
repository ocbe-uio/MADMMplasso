#include <RcppArmadillo.h>
arma::field<arma::cube> admm_MADMMplasso_cpp(
  const arma::vec beta0,
  const arma::mat theta0,
  arma::mat beta,
  arma::mat beta_hat,
  arma::cube theta,
  const double rho1,
  const arma::mat X,
  const arma::mat Z,
  const int max_it,
  const arma::mat W_hat,
  const arma::mat XtY,
  const arma::mat y,
  const int N,
  const double e_abs,
  const double e_rel,
  const double alpha,
  const arma::vec lambda,
  const double alph,
  const arma::mat svd_w_tu,
  const arma::mat svd_w_tv,
  const arma::vec svd_w_d,
  const arma::sp_mat C,
  const arma::vec CW,
  const arma::rowvec gg,
  const bool my_print = true
);
