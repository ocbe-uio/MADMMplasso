#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat model_intercept(const arma::mat beta, const arma::mat X) {
  // The pliable lasso model described in the paper
  // y ~ f(X)
  // formulated as
  // y ~ b_0 + Z theta_0 + X b + \sum( w_j theta_ji )

  arma::mat shared_model = X * beta;
  return shared_model;
}
