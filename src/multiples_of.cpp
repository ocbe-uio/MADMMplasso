#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::ivec multiples_of(arma::ivec x, int divisor) {
  // Returns a vector of booleans indicating whether each element of x is a
  // multiple of divisor.
  for (arma::uword i = 0; i < x.size(); i++) {
    x[i] = (x[i] % divisor == 0)? true: false;
  }
  return x;
}
