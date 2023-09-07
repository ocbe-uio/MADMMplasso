#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec modulo(arma::uvec x, int n) {
    return x - arma::floor(x / n) * n;
}
