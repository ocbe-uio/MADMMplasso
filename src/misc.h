#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>

arma::ivec multiples_of(arma::ivec, int, bool = false);
arma::mat scale_cpp(arma::mat, arma::vec);

#endif
