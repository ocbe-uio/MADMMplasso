# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec sqrt_sum_squared_rows(arma::mat x) {
  // Equivalent to R sqrt(apply(x ^ 2, 1, sum, na.rm = TRUE))
  arma::mat x_sq = arma::pow(x, 2);
  x_sq.replace(arma::datum::nan, 0);
  return arma::sqrt(arma::sum(x_sq, 1));
}
