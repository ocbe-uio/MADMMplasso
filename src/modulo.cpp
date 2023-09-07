#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec modulo(arma::vec x, int n) {
    return x - arma::floor(x / n) * n;
}
