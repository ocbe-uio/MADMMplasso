#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List reg(
  const arma::mat r,
  const arma::mat Z
){

  arma::rowvec beta01(r.n_cols, arma::fill::zeros);
  arma::mat theta01(Z.n_cols, r.n_cols, arma::fill::zeros);

  for (arma::uword e = 0; e < r.n_cols; e++) {
    arma::vec new1 = arma::solve(Z, r.col(e));
    beta01(e) = new1(0);
    theta01.col(e) = new1.tail(new1.n_elem);
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta0") = beta01,
    Rcpp::Named("theta0") = theta01
  );
  return out;
}