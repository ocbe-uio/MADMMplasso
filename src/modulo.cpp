#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::ivec modulo(arma::ivec x, int n) {
    return x - arma::floor(x / n) * n;
}
