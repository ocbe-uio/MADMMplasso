#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int count_nonzero_a_cpp(SEXP x) {
  Rcpp::NumericVector vec(x);
  int n = Rcpp::as<Rcpp::IntegerVector>(vec.attr("dim")).length();

  if (n == 2) {
    arma::mat mat = Rcpp::as<arma::mat>(vec);
    arma::rowvec count1(mat.n_cols, arma::fill::zeros);
    for (unsigned int ww = 0; ww < mat.n_cols; ++ww) {
        count1(ww) = arma::accu(mat.col(ww) != 0);
    }
    return arma::max(count1);
  } else if (n == 3) {
    arma::cube cube = Rcpp::as<arma::cube>(vec);
    arma::vec count1(cube.n_slices, arma::fill::zeros);
    for (unsigned int ww = 0; ww < cube.n_slices; ++ww) {
        count1(ww) = arma::accu(cube.slice(ww) != 0);
    }
    return arma::max(count1);
  } else {
    Rcpp::stop("Input must be a 2D or 3D matrix.");
  }
}
