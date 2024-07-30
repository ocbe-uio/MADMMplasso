#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec lm_arma(const arma::vec &R, const arma::mat &Z) {
    // Add a column of ones to Z
    arma::mat Z_intercept = arma::join_rows(arma::ones<arma::vec>(Z.n_rows), Z);

    // Solve the system of linear equations
    arma::vec coefficients = arma::solve(Z_intercept, R);

    // Replace 0 with NA (arma::datum::nan)
    for (arma::uword i = 0; i < coefficients.n_elem; i++) {
        if (coefficients(i) == 0) {
            coefficients(i) = arma::datum::nan;
        }
    }

    return coefficients;
}

// [[Rcpp::export]]
arma::mat reg(const arma::mat r, const arma::mat Z) {
  arma::rowvec beta01(r.n_cols, arma::fill::zeros);
  arma::mat theta01(Z.n_cols, r.n_cols, arma::fill::zeros);

  for (arma::uword e = 0; e < r.n_cols; e++) {
    arma::vec new1 = lm_arma(r.col(e), Z);
    beta01(e) = new1(0);
    theta01.col(e) = new1.tail(new1.n_elem - 1);
  }

  arma::mat out = arma::join_vert(beta01, theta01);
  return out;
}
