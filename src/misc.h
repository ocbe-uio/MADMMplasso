#include <RcppArmadillo.h>
arma::ivec multiples_of(arma::ivec, int, bool = false);
arma::mat scale_cpp(arma::mat, arma::vec);
arma::vec sqrt_sum_squared_rows(arma::mat);
arma::uvec modulo(arma::uvec, int);
arma::mat model_intercept(
  const arma::vec,
  const arma::mat,
  const arma::mat,
  const arma::cube,
  const arma::mat,
  const arma::mat
);
arma::mat model_p(
  const arma::vec,
  const arma::mat,
  const arma::mat,
  const arma::cube,
  const arma::mat,
  const arma::mat
);
arma::mat reg(const arma::mat, const arma::mat);
int count_nonzero_a_cpp(SEXP);
int count_nonzero_a_sp_mat(arma::sp_mat);
int count_nonzero_a_cube(arma::cube);
int count_nonzero_a_mat(arma::mat);
