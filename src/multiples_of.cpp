#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::ivec multiples_of(arma::ivec x, const int divisor, const bool subset_out = false) {
  // Returns a vector of the elements of x that are a multiple of divisor.

  // Validation
  if (divisor == 0) {
    Rcpp::stop("divisor cannot be 0");
  }

  // Finding and returning multiples
  arma::ivec div_idx = arma::zeros<arma::ivec>(x.size());
  for (arma::uword i = 0; i < x.size(); i++) {
    if (x[i] % divisor == 0) {
      div_idx[i] = i + 1;
    } else {
      div_idx[i] = 0;
      x[i] = 0;
    }
  }

  // Returning either nonzero x or div_idx depending on subset_out
  arma::ivec out = (subset_out) ? arma::nonzeros(x) : arma::nonzeros(div_idx);
  return out;
}
