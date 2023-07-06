// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// admm_MADMMplasso_cpp
Rcpp::List admm_MADMMplasso_cpp(arma::vec beta0, arma::mat theta0, arma::mat beta, arma::mat beta_hat, arma::cube theta, double rho1, arma::mat X, arma::mat Z, int max_it, arma::mat W_hat, arma::mat XtY, arma::mat y, int N, int p, int K, double e_abs, double e_rel, double alpha, arma::vec lambda, double alph, Rcpp::List svd_w, Rcpp::List tree, Rcpp::List invmat, arma::cube V, arma::cube Q, arma::mat E, arma::cube EE, arma::cube O, arma::cube P, arma::mat H, arma::cube HH, arma::vec gg, bool my_print);
RcppExport SEXP _MADMMplasso_admm_MADMMplasso_cpp(SEXP beta0SEXP, SEXP theta0SEXP, SEXP betaSEXP, SEXP beta_hatSEXP, SEXP thetaSEXP, SEXP rho1SEXP, SEXP XSEXP, SEXP ZSEXP, SEXP max_itSEXP, SEXP W_hatSEXP, SEXP XtYSEXP, SEXP ySEXP, SEXP NSEXP, SEXP pSEXP, SEXP KSEXP, SEXP e_absSEXP, SEXP e_relSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP alphSEXP, SEXP svd_wSEXP, SEXP treeSEXP, SEXP invmatSEXP, SEXP VSEXP, SEXP QSEXP, SEXP ESEXP, SEXP EESEXP, SEXP OSEXP, SEXP PSEXP, SEXP HSEXP, SEXP HHSEXP, SEXP ggSEXP, SEXP my_printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W_hat(W_hatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtY(XtYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type e_abs(e_absSEXP);
    Rcpp::traits::input_parameter< double >::type e_rel(e_relSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type svd_w(svd_wSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type invmat(invmatSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type E(ESEXP);
    Rcpp::traits::input_parameter< arma::cube >::type EE(EESEXP);
    Rcpp::traits::input_parameter< arma::cube >::type O(OSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type HH(HHSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gg(ggSEXP);
    Rcpp::traits::input_parameter< bool >::type my_print(my_printSEXP);
    rcpp_result_gen = Rcpp::wrap(admm_MADMMplasso_cpp(beta0, theta0, beta, beta_hat, theta, rho1, X, Z, max_it, W_hat, XtY, y, N, p, K, e_abs, e_rel, alpha, lambda, alph, svd_w, tree, invmat, V, Q, E, EE, O, P, H, HH, gg, my_print));
    return rcpp_result_gen;
END_RCPP
}
// multiples_of
arma::ivec multiples_of(arma::ivec x, int divisor, bool subset_out);
RcppExport SEXP _MADMMplasso_multiples_of(SEXP xSEXP, SEXP divisorSEXP, SEXP subset_outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< bool >::type subset_out(subset_outSEXP);
    rcpp_result_gen = Rcpp::wrap(multiples_of(x, divisor, subset_out));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MADMMplasso_admm_MADMMplasso_cpp", (DL_FUNC) &_MADMMplasso_admm_MADMMplasso_cpp, 33},
    {"_MADMMplasso_multiples_of", (DL_FUNC) &_MADMMplasso_multiples_of, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MADMMplasso(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
