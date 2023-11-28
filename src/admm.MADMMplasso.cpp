#include <RcppArmadillo.h>
#include "misc.h"

// [[Rcpp::depends(RcppArmadillo)]]
// #' Fit the ADMM part of  model for a given lambda vale
// #' @param X  n by p matrix of predictors
// #' @param Z n by nz matrix of modifying variables. The elements of z
// #' may represent quantitative or categorical variables, or a mixture of the two.
// #' Categorical varables should be coded by 0-1 dummy variables: for a k-level
// #' variable, one can use either k or k-1  dummy variables.
// #' @param beta0 TODO: fill paramater description
// #' @param theta0 TODO: fill paramater description
// #' @param beta TODO: fill paramater description
// #' @param beta_hat TODO: fill paramater description
// #' @param theta TODO: fill paramater description
// #' @param rho1 TODO: fill paramater description
// #' @param max_it TODO: fill paramater description
// #' @param W_hat TODO: fill paramater description
// #' @param XtY TODO: fill paramater description
// #' @param y TODO: fill paramater description
// #' @param N TODO: fill paramater description
// #' @param p TODO: fill paramater description
// #' @param K TODO: fill paramater description
// #' @param e.abs TODO: fill paramater description
// #' @param e.rel TODO: fill paramater description
// #' @param alpha TODO: fill paramater description
// #' @param lambda TODO: fill paramater description
// #' @param alph TODO: fill paramater description
// #' @param svd.w TODO: fill paramater description
// #' @param tree TODO: fill paramater description
// #' @param my_print TODO: fill paramater description
// #' @param invmat TODO: fill paramater description
// #' @param V TODO: fill paramater description
// #' @param Q TODO: fill paramater description
// #' @param E TODO: fill paramater description
// #' @param EE TODO: fill paramater description
// #' @param O TODO: fill paramater description
// #' @param P TODO: fill paramater description
// #' @param H TODO: fill paramater description
// #' @param HH TODO: fill paramater description
// #' @param gg TODO: fill paramater description
// #' @return  predicted values for the ADMM part
// #' @description TODO: add description
// #' @export
// [[Rcpp::export]]
Rcpp::List admm_MADMMplasso_cpp(
  const arma::vec beta0,
  const arma::mat theta0,
  arma::mat beta,
  arma::mat beta_hat,
  arma::cube theta,
  const double rho1,
  const arma::mat X,
  const arma::mat Z,
  const int max_it,
  const arma::mat W_hat,
  const arma::mat XtY,
  const arma::mat y,
  const int N,
  const double e_abs,
  const double e_rel,
  const double alpha,
  const arma::vec lambda,
  const double alph,
  const Rcpp::List svd_w,
  const Rcpp::List tree,
  const Rcpp::List invmat,
  const arma::vec gg,
  const bool my_print = true
) {
  const Rcpp::List TT = tree;
  const arma::sp_mat C = TT["Tree"];
  const arma::vec CW = TT["Tw"];
  const arma::mat svd_w_tu = Rcpp::as<arma::mat>(svd_w["u"]).t();
  const arma::mat svd_w_tv = Rcpp::as<arma::mat>(svd_w["v"]).t();
  const int D = y.n_cols;
  const int p = X.n_cols;
  const uint K = Z.n_cols;

  arma::cube V(p, 2 * (1 + K), D, arma::fill::zeros);
  arma::cube O(p, 2 * (1 + K), D, arma::fill::zeros);
  arma::mat E(y.n_cols * C.n_rows, p + p * K);
  arma::cube EE(p, 1 + K, D, arma::fill::zeros);

  arma::cube Q(p, 1 + K, D, arma::fill::zeros);
  arma::cube P(p, 1 + K, D, arma::fill::zeros);
  arma::mat H(y.n_cols * C.n_rows, p + p * K);
  arma::cube HH(p, 1 + K, D, arma::fill::zeros);


  // for response groups =======================================================
  const arma::ivec input = Rcpp::seq_len(D * C.n_rows);
  arma::mat I = arma::zeros<arma::mat>(C.n_rows * D, D);
  arma::ivec II = multiples_of(input, D, true);

  // Updating I ================================================================
  I.rows(0, D - 1).diag() = arma::conv_to<arma::mat>::from(C.row(0) * CW(0));

  int c_count = 2;
  for (unsigned int e = 0; e < II.n_elem - 1; ++e) {
    I.rows(II(e), c_count * D - 1).diag() = arma::conv_to<arma::mat>::from(C.row(c_count - 1) * CW(c_count - 1));
    ++c_count;
  }
  const arma::vec new_I = diagvec(I.t() * I);

  // Overlapping group from covariates =========================================
  arma::umat G = arma::zeros<arma::umat>(2 * (1 + K), 1 + K);
  const arma::umat diag_G = arma::eye<arma::umat>(K + 1, K + 1);
  for (uint i = 0; i < K + 1; ++i) {
    G.row(i) = diag_G.row(i);
  }
  for (uint i = K + 2; i < 2 * (K + 1); ++i) {
    G.row(i) = diag_G.row(i - (K + 1));
  }

  arma::cube V_old = V;
  arma::cube Q_old = Q;
  arma::mat E_old = E;
  arma::cube EE_old = EE;
  double res_pri = 0.;
  double res_dual = 0.;
  const arma::mat SVD_D = arma::diagmat(Rcpp::as<arma::vec>(svd_w["d"]));
  const arma::mat R_svd = (svd_w_tu.t() * SVD_D) / N;
  double rho = rho1;
  arma::cube Big_beta11 = V;
  arma::mat res_val;  // declared here because it's also needed outside the loop
  bool converge = false;
  if (my_print) {
    Rcpp::Rcout << "\ni\tres_dual\te_dual\t\tres_pri\t\te_primal" << std::endl;
  }
  for (int i = 1; i < max_it + 1; i++) {
    arma::mat shared_model = model_intercept(beta0, theta0, beta_hat, theta, W_hat, Z);
    arma::mat r_current = y - shared_model;
    Rcpp::List b = reg(r_current, Z);
    arma::mat beta0 = b["beta0"];
    arma::mat theta0 = b["theta0"];
    arma::mat new_y = y - (arma::ones(N) * beta0 + Z * theta0);
    arma::mat XtY = W_hat.t() * new_y;
    arma::cube main_beta(p, K + 1, D, arma::fill::zeros);
    res_val = rho * (I.t() * E - (I.t() * H));
    arma::vec v_diff1(D, arma::fill::zeros);
    arma::vec q_diff1(D, arma::fill::zeros);
    arma::vec ee_diff1(D, arma::fill::zeros);

    arma::vec new_G(p + p * K, arma::fill::zeros);
    new_G.rows(0, p - 1).fill(1);
    new_G.rows(p, p + p * K - 1).fill(2);
    new_G = rho * (1 + new_G);
    arma::cube invmat(new_G.n_rows, 1, D);  // denominator of the beta estimates
    for (arma::uword slc = 0; slc < D; slc++) {
      invmat.slice(slc) = new_G + rho * (new_I(slc) + 1);
    }

    for (int rr = 0; rr < D; rr++) {
      double DD1 = rho * (new_I(rr) + 1);
      arma::vec DD2 = new_G + DD1;
      invmat.slice(rr) = DD2;  // Matrix::chol2inv( Matrix::chol(new_sparse) )
    }

    for (int jj = 0; jj < D; jj++) {
      arma::mat group = rho * (G.t() * V.slice(jj).t() - G.t() * O.slice(jj).t());
      arma::vec group1 = group.row(0).t();
      arma::mat group2 = group.tail_rows(group.n_rows - 1).t();
      arma::mat new_group(p, K + 1, arma::fill::zeros);
      new_group.col(0) = group1;
      new_group.tail_cols(new_group.n_cols - 1) = group2;
      arma::vec my_beta_jj = XtY.col(jj) / N +\
        arma::vectorise(new_group) + res_val.row(jj).t() +\
        arma::vectorise(rho * (Q.slice(jj) - P.slice(jj))) +\
        arma::vectorise(rho * (EE.slice(jj) - HH.slice(jj)));
      arma::mat DD3 = arma::diagmat(1 / invmat.slice(jj));
      arma::mat part_z = DD3 * W_hat.t();
      arma::mat part_y = DD3 * my_beta_jj;

      arma::mat beta_hat_j = arma::inv(arma::inv(R_svd) + svd_w_tv * part_z);
      beta_hat_j = beta_hat_j * (svd_w_tv * part_y);
      beta_hat_j = part_z * beta_hat_j;

      arma::vec beta_hat_JJ = arma::vectorise(part_y - beta_hat_j);
      beta_hat.col(jj) = beta_hat_JJ;
      arma::mat beta_hat1 = arma::reshape(beta_hat_JJ, p, 1 + K);
      arma::mat b_hat = alph * beta_hat1 + (1 - alph) * Q.slice(jj);
      Q.slice(jj).col(0) = b_hat.col(0) + P.slice(jj).col(0);
      arma::mat new_mat = b_hat.tail_cols(b_hat.n_cols - 1) + P.slice(jj).tail_cols(P.slice(jj).n_cols - 1);
      Q.slice(jj).tail_cols(Q.n_cols - 1) = arma::sign(new_mat) % arma::max(arma::abs(new_mat) - ((alpha * lambda(jj)) / rho), arma::zeros(arma::size(new_mat)));
      b_hat = alph * beta_hat1 + (1 - alph) * EE.slice(jj);
      new_mat = b_hat + HH.slice(jj);

      arma::vec row_norm1 = sqrt_sum_squared_rows(new_mat);

      arma::vec coef_term1 = arma::max(1 - gg(1) / rho / row_norm1, arma::zeros(arma::size(row_norm1)));

      arma::mat ee1 = scale_cpp(new_mat.t(), 1 / coef_term1);
      EE.slice(jj) = arma::trans(ee1);

      arma::mat Big_beta = arma::trans(G * beta_hat1.t());
      Big_beta11.slice(jj) = Big_beta;
      arma::mat Big_beta1 = alph * Big_beta + (1 - alph) * V.slice(jj);

      // Now we have the main part.
      new_mat = Big_beta1 + O.slice(jj);
      arma::mat new_mat1 = new_mat.head_cols(K + 1);
      arma::mat new_mat2 = new_mat.tail_cols(new_mat.n_cols - K - 1);
      row_norm1 = sqrt_sum_squared_rows(new_mat1);
      arma::vec row_norm2 = sqrt_sum_squared_rows(new_mat2);
      coef_term1 = arma::max(1 - (1 - alpha) * lambda(jj) / rho / row_norm1, arma::zeros(arma::size(row_norm1)));
      arma::vec coef_term2 = arma::max(1 - (1 - alpha) * lambda(jj) / rho / row_norm2, arma::zeros(arma::size(row_norm2)));
      arma::mat N_V1 = scale_cpp(new_mat1.t(), 1 / coef_term1);
      arma::mat N_V2 = scale_cpp(new_mat2.t(), 1 / coef_term2);
      V.slice(jj) = arma::join_horiz(N_V1.t(), N_V2.t());
      P.slice(jj) += beta_hat1 - Q.slice(jj);
      HH.slice(jj) += beta_hat1 - EE.slice(jj);
      O.slice(jj) += Big_beta - V.slice(jj);
      v_diff1(jj) = arma::accu(arma::pow(Big_beta - V.slice(jj), 2));
      q_diff1(jj) = arma::accu(arma::pow(beta_hat1 - Q.slice(jj), 2));
      ee_diff1(jj) = arma::accu(arma::pow(beta_hat1 - EE.slice(jj), 2));
    }

    arma::mat Big_beta_response = I * beta_hat.t();
    arma::mat b_hat_response = alph * Big_beta_response + (1 - alph) * E;
    arma::mat new_mat = b_hat_response + H;
    arma::cube new_mat_group(p + p * K, y.n_cols, C.n_rows);
    arma::cube beta_group(p + p * K, y.n_cols, C.n_rows);
    new_mat_group.slice(0) = new_mat.rows(0, y.n_cols - 1).t();
    beta_group.slice(0) = Big_beta_response.rows(0, y.n_cols - 1).t();
    arma::mat beta_transform(p, (K + 1) * y.n_cols, arma::fill::zeros);
    beta_transform.cols(0, K) = arma::reshape(new_mat_group.slice(0).col(0), p, K + 1);
    arma::uvec input2 = arma::regspace<arma::uvec>(1, y.n_cols * (1 + K));
    arma::uvec multiple_of_K = modulo(input2, K + 1);
    arma::uvec II2 = input2.elem(arma::find(multiple_of_K == 0));
    int e2 = II2(0);

    for (arma::uword c_count2 = 2; c_count2 < y.n_cols + 1; c_count2++) {
      arma::mat new_mat_group_sub = arma::trans(new_mat_group.slice(0).col(c_count2 - 1));
      beta_transform.cols(e2, (c_count2 * (1 + K)) - 1) = new_mat_group_sub.reshape(p, K + 1);
      e2 = II2(c_count2 - 1);
    }

    arma::vec norm_res = sqrt_sum_squared_rows(beta_transform);
    arma::vec coef_term1 = arma::max(1 - gg(0) / rho / norm_res, arma::zeros(arma::size(norm_res)));
    arma::mat N_E1 = scale_cpp(beta_transform.t(), 1 / coef_term1).t();
    arma::mat beta_transform1(p + p * K, y.n_cols, arma::fill::zeros);
    beta_transform1.col(0) = N_E1.head_cols(K + 1).as_col();
    arma::uvec input3 = arma::regspace<arma::uvec>(1, y.n_cols * (1 + K));
    multiple_of_K = modulo(input3, K + 1);
    arma::uvec II3 = input2.elem(arma::find(multiple_of_K == 0));
    int e3 = II3(0);

    for (arma::uword c_count3 = 1; c_count3 < y.n_cols; c_count3++) {
      beta_transform1.col(c_count3) = N_E1.cols(e3, (K + 1) * (c_count3 + 1) - 1).as_col();
      e3 = II3(c_count3);
    }

    // Original code has N_E as list, but it's only storing beta_transform1.t()
    arma::cube N_E(beta_transform1.n_cols, beta_transform1.n_rows, C.n_rows);

    N_E.slice(0) = beta_transform1.t();
    int e = II(0);

    for (arma::uword c_count = 1; c_count < C.n_rows; c_count++) {
      new_mat_group.slice(c_count) = new_mat.rows(e, ((c_count + 1) * y.n_cols) - 1).t();
      beta_group.slice(c_count) = Big_beta_response.rows(e, ((c_count + 1) * y.n_cols) - 1).t();
      arma::mat beta_transform(p, (K + 1) * y.n_cols, arma::fill::zeros);
      beta_transform.cols(0, K) = arma::reshape(new_mat_group.slice(c_count).col(0), p, K + 1);
      arma::uvec input2 = arma::regspace<arma::uvec>(1, y.n_cols * (1 + K));
      arma::uvec multiple_of_K = modulo(input2, K + 1);
      arma::uvec II2 = input2.elem(arma::find(multiple_of_K == 0));
      int e2 = II2(0);

      for (arma::uword c_count2 = 1; c_count2 < y.n_cols; c_count2++) {
        beta_transform.cols(e2, c_count2 * (K + 1) + K) = arma::reshape(new_mat_group.slice(c_count).col(c_count2), p, K + 1);
        e2 = II2(c_count2);
      }

      norm_res = sqrt_sum_squared_rows(beta_transform);
      coef_term1 = arma::max(1 - gg(0) / rho / norm_res, arma::zeros(arma::size(norm_res)));
      N_E1 = scale_cpp(beta_transform.t(), 1 / coef_term1).t();
      beta_transform1 = arma::zeros<arma::mat>(p + p * K, y.n_cols);
      beta_transform1.col(0) = N_E1.head_cols(K + 1).as_col();
      arma::uvec input3 = arma::regspace<arma::uvec>(1, y.n_cols * (1 + K));
      multiple_of_K = modulo(input3, K + 1);
      arma::uvec II3 = input2.elem(arma::find(multiple_of_K == 0));
      int e3 = II3(0);

      for (arma::uword c_count3 = 1; c_count3 < y.n_cols; c_count3++) {
        beta_transform1.col(c_count3) = N_E1.cols(e3, (K + 1) * (c_count3 + 1) - 1).as_col();
        e3 = II3(c_count3);
      }

      N_E.slice(c_count) = beta_transform1.t();
      e = II(c_count);
    }

      E.rows(0, C.n_cols - 1) = N_E.slice(0);
      e = II(0);

      for (arma::uword c_count = 1; c_count < C.n_rows; c_count++) {
        E.rows(e, ((c_count + 1) * y.n_cols) - 1) = N_E.slice(c_count);
        e = II(c_count);
      }

      H += Big_beta_response - E;

      double v_diff = arma::accu(arma::pow(-rho * (V - V_old), 2));
      double q_diff = arma::accu(arma::pow(-rho * (Q - Q_old), 2));
      double e_diff = arma::accu(arma::pow(-rho * (E - E_old), 2));
      double ee_diff = arma::accu(arma::pow(-rho * (EE - EE_old), 2));
      double s = sqrt(v_diff + q_diff + e_diff + ee_diff);

      double v_diff1_sum = arma::accu(v_diff1);
      double q_diff1_sum = arma::accu(q_diff1);
      double e_diff1 = arma::accu(arma::pow(Big_beta_response - E, 2));
      double ee_diff1_sum = arma::accu(ee_diff1);
      double r = sqrt(v_diff1_sum + q_diff1_sum + e_diff1 + ee_diff1_sum);

      res_dual = s;
      res_pri = r;

      double part_1 = sqrt(Big_beta11.n_elem + 2 * beta_hat.n_elem + Big_beta_response.n_elem);
      arma::vec Big_beta11_vec = arma::vectorise(Big_beta11);
      arma::vec beta_hat_vec = arma::vectorise(beta_hat);
      arma::vec Big_beta_response_vec = arma::vectorise(Big_beta_response);
      arma::vec part_2 = arma::join_vert(Big_beta11_vec, beta_hat_vec, beta_hat_vec, Big_beta_response_vec);
      double part_2_norm = arma::norm(part_2);

      arma::vec V_vec = arma::vectorise(V);
      arma::vec Q_vec = arma::vectorise(Q);
      arma::vec E_vec = arma::vectorise(E);
      arma::vec EE_vec = arma::vectorise(EE);
      arma::vec part_3 = -1 * arma::join_vert(V_vec, Q_vec, E_vec, EE_vec);
      double part_3_norm = arma::norm(part_3);

      double part_2_3 = std::max(part_2_norm, part_3_norm);

      double e_primal = part_1 * e_abs + e_rel * part_2_3;

      arma::vec O_vec = arma::vectorise(O);
      arma::vec P_vec = arma::vectorise(P);
      arma::vec H_vec = arma::vectorise(H);
      arma::vec HH_vec = arma::vectorise(HH);
      arma::vec part_4 = arma::join_vert(O_vec, P_vec, H_vec, HH_vec);
      double e_dual = part_1 * e_abs + e_rel * arma::max(part_4);

      V_old = V,
      Q_old = Q,
      E_old = E,
      EE_old = EE;

      if (res_pri > 10 * res_dual) {
        rho *= 2;
      } else if (res_pri * 10 < res_dual ) {
        rho /= 2;
      }

      if (my_print) {
        printf("%u\t%e\t%e\t%e\t%e\n", i, res_dual, e_dual, res_pri, e_primal);
      }

      if (res_pri <= e_primal && res_dual <= e_dual) {
        if (my_print) {
          printf("Convergence reached after %u iterations", i);
        }
        converge = true;
        break;
      }
  }

  res_val = I.t() * E;

  for (arma::uword jj = 0; jj < y.n_cols; jj++) {
    arma::mat group = G.t() * V.slice(jj).t();
    arma::vec group1 = group.row(0).t();
    arma::mat group2 = group.tail_rows(group.n_rows - 1).t();
    arma::mat new_group = arma::zeros<arma::mat>(p, K + 1);
    new_group.col(0) = group1;
    new_group.tail_cols(new_group.n_cols - 1) = group2;
    arma::vec new_g_theta = arma::vectorise(new_group);

    arma::mat beta_hat_B1 = beta_hat.submat(0, jj, p - 1, jj);
    arma::vec new_g_theta_B1 = new_g_theta.subvec(0, p - 1);
    arma::vec Q_B1 = Q.slice(jj).col(0);
    arma::vec finB1 = arma::vectorise(beta_hat_B1 % (new_g_theta_B1 != 0) % (Q_B1 != 0));

    arma::mat beta_hat_B2 = beta_hat.submat(p, jj, beta_hat.n_rows - 1, jj);
    arma::vec new_g_theta_B2 = new_g_theta.subvec(p, new_g_theta.n_elem - 1);
    arma::vec Q_B2 = arma::vectorise(Q.subcube(0, 1, jj, Q.n_rows - 1, Q.n_cols - 1, jj));
    arma::vec finB2 = arma::vectorise(beta_hat_B2 % (new_g_theta_B2 != 0) % (Q_B2 != 0));

    arma::mat beta_hat1 = arma::reshape(arma::join_vert(finB1, finB2), p, K + 1);
    beta.col(jj) = beta_hat1.col(0);
    theta.slice(jj) = beta_hat1.tail_cols(beta_hat1.n_cols - 1);
    beta_hat.col(jj) = arma::join_vert(beta_hat1.col(0), arma::vectorise(theta.slice(jj)));
  }

  arma::mat y_hat = model_p(beta0, theta0, beta_hat, theta, W_hat, Z);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("theta0") = theta0,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("converge") = converge,
    Rcpp::Named("obj") = NULL,
    Rcpp::Named("beta_hat") = beta_hat,
    Rcpp::Named("y_hat") = y_hat
  );
  out.attr("class") = "admm.MADMMplasso";
  return out;
}
