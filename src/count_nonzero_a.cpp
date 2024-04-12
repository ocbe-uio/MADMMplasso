#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int count_nonzero_a_cpp(SEXP x) { // FIXME: doesn't work with arma types
  // This function counts the maximum count of of non-zero elements in
  // the columns of a matrix or the slices of a cube
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

// [[Rcpp::export]]
int count_nonzero_a_sp_mat(arma::sp_mat x) {
  // This function counts the maximum count of of non-zero elements in
  // the columns of a matrix or the slices of a cube
  arma::mat x_dense = arma::mat(x);
  arma::rowvec count1(x_dense.n_cols, arma::fill::zeros);
  for (unsigned int ww = 0; ww < x_dense.n_cols; ++ww) {
    count1(ww) = arma::accu(x_dense.col(ww) != 0);
  }
  return arma::max(count1);
}

// [[Rcpp::export]]
int count_nonzero_a_cube(arma::cube x) {
  arma::vec count1(x.n_slices, arma::fill::zeros);
  for (unsigned int ww = 0; ww < x.n_slices; ++ww) {
    count1(ww) = arma::accu(x.slice(ww) != 0);
  }
  return arma::max(count1);
}

// [[Rcpp::export]]
int count_nonzero_a_mat(arma::mat x) {
  arma::vec count1(x.n_cols, arma::fill::zeros);
  for (unsigned int ww = 0; ww < x.n_cols; ++ww) {
    count1(ww) = arma::accu(x.col(ww) != 0);
  }
  return arma::max(count1);
}
