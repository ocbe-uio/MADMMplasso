#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat model_p(
  const arma::vec beta0,
  const arma::mat theta0,
  const arma::mat beta,
  const arma::mat X,
  const arma::mat Z
){
  arma::mat intercepts = arma::ones(X.n_rows) * beta0.t() + Z * theta0;
  arma::mat shared_model = X * beta;
  return(intercepts + shared_model);
}
