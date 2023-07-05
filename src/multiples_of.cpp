#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::ivec multiples_of(arma::ivec x, int divisor) {
  // Returns a vector of the elements of x that are a multiple of divisor.

  // Validation
  if (divisor == 0) {
    Rcpp::stop("divisor cannot be 0");
  }

  // Finding and returning multiples
  arma::ivec div_idx = arma::zeros<arma::ivec>(x.size());
  for (arma::uword i = 0; i < x.size(); i++) {
    div_idx[i] = (x[i] % divisor == 0)? i + 1: 0;
  }
  return arma::nonzeros(div_idx);
}
