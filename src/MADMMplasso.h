#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>
Rcpp::List admm_MADMMplasso_cpp(
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
  const Rcpp::List svd_w,
  const Rcpp::List tree,
  const Rcpp::List invmat,
  const arma::vec gg,
  const bool my_print = true
);
#endif
