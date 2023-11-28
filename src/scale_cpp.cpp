#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat scale_cpp(arma::mat x, arma::vec scale) {
  // ad-hoc implementation of R's scale function. Works if center = FALSE and
  // scale is a vector with length equal to the number of columns of x.
  return arma::trans(x.t() / arma::repmat(scale, 1, x.n_rows));
}
