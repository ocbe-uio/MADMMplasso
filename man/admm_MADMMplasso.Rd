% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/admm_MADMMplasso.R
\name{admm_MADMMplasso}
\alias{admm_MADMMplasso}
\title{Fit the ADMM part of  model for the given lambda values}
\usage{
admm_MADMMplasso(
  beta0,
  theta0,
  beta,
  beta_hat,
  theta,
  rho1,
  X,
  Z,
  max_it,
  W_hat,
  XtY,
  y,
  N,
  e.abs,
  e.rel,
  alpha,
  lambda,
  alph,
  svd.w,
  tree,
  my_print,
  invmat,
  gg = 0.2
)
}
\arguments{
\item{beta0}{a vector of length ncol(y) of estimated beta_0 coefficients}

\item{theta0}{matrix of the initial theta_0 coefficients  ncol(Z) by ncol(y)}

\item{beta}{a matrix of the initial beta coefficients    ncol(X) by ncol(y)}

\item{beta_hat}{a matrix of the initial beta and theta coefficients    (ncol(X)+ncol(X) by ncol(Z)) by ncol(y)}

\item{theta}{an array of initial theta coefficients    ncol(X) by ncol(Z) by ncol(y)}

\item{rho1}{the Lagrange variable for the ADMM which is usually included as rho in the MADMMplasso call.}

\item{X}{N by p matrix of predictors}

\item{Z}{N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
Categorical variables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.}

\item{max_it}{maximum number of iterations in loop for one lambda during the ADMM optimization}

\item{W_hat}{N by (p+(p by nz)) of the main and interaction predictors. This generated internally  when MADMMplasso is called or by using the function generate_my_w.}

\item{XtY}{a matrix formed by multiplying the transpose of X by y.}

\item{y}{N by D matrix  of responses. The X and Z variables are centered in the function. We recommend that X and Z also be standardized before the call}

\item{N}{nrow(X)}

\item{e.abs}{absolute error for the ADMM}

\item{e.rel}{relative error for the ADMM}

\item{alpha}{mixing parameter. When the goal is to include more interactions, alpha should be very small and vice versa.}

\item{lambda}{user specified lambda_3 values.}

\item{alph}{an overrelaxation parameter in [1, 1.8]. The implementation is borrowed from Stephen Boyd's \href{https://stanford.edu/~boyd/papers/admm/lasso/lasso.html}{MATLAB code}}

\item{svd.w}{singular value decomposition of W}

\item{tree}{The results from the hierarchical clustering of the response matrix. The easy way to obtain this is by using the function (tree_parms) which gives a default clustering. However, user decide on a specific structure and then input a tree that follows such structure.}

\item{my_print}{Should information form each ADMM iteration be printed along the way? This prints the dual and primal residuals}

\item{invmat}{A list of length ncol(y), each containing the C_d part of equation 32 in the paper}

\item{gg}{penalty terms for the tree structure for lambda_1 and  lambda_2 for the admm call.}
}
\value{
predicted values for the ADMM part
beta0:  estimated beta_0 coefficients  having a size of 1 by ncol(y)

beta: estimated beta coefficients  having a matrix   ncol(X) by ncol(y)

BETA_hat:  estimated beta and theta coefficients  having a matrix   (ncol(X)+ncol(X) by ncol(Z)) by ncol(y)

theta0:  estimated theta_0 coefficients  having a matrix   ncol(Z) by ncol(y)

theta:  estimated theta coefficients  having a an array   ncol(X) by ncol(Z) by ncol(y)
converge: did the algorithm converge?

Y_HAT:  predicted response nrow(X) by ncol(y)
}
\description{
This function fits a multi-response pliable lasso model over a path of regularization values.
}
