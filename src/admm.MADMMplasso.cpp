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
  const arma::mat beta,
  arma::mat beta_hat,
  const arma::cube theta,
  const double rho1,
  const arma::mat X,
  const arma::mat Z,
  const int max_it,
  const arma::mat W_hat,
  const arma::mat XtY,
  const arma::mat y,
  const int N,
  const int p,
  const uint K,
  const double e_abs,
  const double e_rel,
  const double alpha,
  const arma::vec lambda,
  const double alph,
  const Rcpp::List svd_w,
  const Rcpp::List tree,
  const Rcpp::List invmat,
  arma::cube V,
  arma::cube Q,
  const arma::mat E,
  arma::cube EE,
  arma::cube O,
  arma::cube P,
  const arma::mat H,
  arma::cube HH,
  const arma::vec gg,
  const bool my_print = true
) {
  const Rcpp::List TT = tree;
  const arma::sp_mat C = TT["Tree"];
  const arma::vec CW = TT["Tw"];
  const arma::mat svd_w_tu = Rcpp::as<arma::mat>(svd_w["u"]).t();
  const arma::mat svd_w_tv = Rcpp::as<arma::mat>(svd_w["v"]).t();
  const int D = y.n_cols;

  // for response groups =======================================================
  const arma::ivec input = Rcpp::seq_len(D * C.n_rows);
  arma::mat I = arma::zeros<arma::mat>(C.n_rows * D, D);
  const arma::ivec II = multiples_of(input, D, true);

  // Updating I ================================================================
  I.rows(0, D - 1).diag() = arma::conv_to<arma::mat>::from(C.row(0) * CW(0));

  int c_count = 2;
  for (unsigned int e = 0; e < II.n_elem - 1; ++e) {
    I.rows(II(e), c_count * D - 1).diag() = arma::conv_to<arma::mat>::from(C.row(c_count - 1) * CW(c_count - 1));
    ++c_count;
  }
  const arma::vec new_I = diagvec(I.t() * I);

  // Overlapping group from covariates =========================================
  arma::mat G = arma::zeros<arma::mat>(2 * (1 + K), 1 + K);
  const arma::mat diag_G = arma::eye(K + 1, K + 1);
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
  const arma::vec obj;
  const arma::mat SVD_D = arma::diagmat(Rcpp::as<arma::vec>(svd_w["d"]));
  const arma::mat R_svd = (svd_w_tu * SVD_D) / N;
  double rho = rho1;
  arma::cube Big_beta11 = V;

  // Importing R functions (this adds compute overhead)
  // Ideally, these functions should also be ported to C++ to reduce
  // cross-language communication
  Rcpp::Environment MADMMplasso = Rcpp::Environment::namespace_env("MADMMplasso");
  Rcpp::Function model_intercept = MADMMplasso["model_intercept"];
  Rcpp::Function reg = MADMMplasso["reg"];
  Rcpp::Function model_p = MADMMplasso["model_p"];

  bool converge = false;
  for (int i = 1; i < max_it; i++) {
    arma::mat shared_model = Rcpp::as<arma::mat>(model_intercept(beta0, theta0, beta_hat, theta, W_hat, Z));
    arma::mat r_current = y - shared_model;
    Rcpp::List b = reg(r_current,Z);
    arma::mat beta0 = b["beta0"];
    arma::mat theta0 = b["theta0"];
    arma::mat new_y = y - (arma::ones(N) * beta0 + Z * theta0);
    arma::mat XtY = W_hat.t() * new_y;
    arma::cube main_beta(p, K + 1, D, arma::fill::zeros);
    arma::mat res_val = rho * (I.t() * E - (I.t() * H));
    arma::vec v_diff1(D, arma::fill::zeros);
    arma::vec q_diff1(D, arma::fill::zeros);
    arma::vec ee_diff1(D, arma::fill::zeros);
    arma::vec new_G(p + p * K, arma::fill::zeros);
    new_G.rows(0, p - 1).fill(1);
    new_G.rows(p, p + p * K - 1).fill(2);
    new_G = rho * (1 + new_G);
    arma::cube invmat(new_G.n_rows, 1, D);  // denominator of the beta estimates

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
      Q.slice(jj).tail_cols(Q.slice(jj).n_cols - 1) = arma::sign(new_mat) % arma::max(arma::abs(new_mat) - ((alpha * lambda(jj)) / rho), arma::zeros(arma::size(new_mat)));
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
    Rcpp::List N_E;
    new_mat_group.slice(0) = new_mat.rows(0, y.n_cols - 1).t();
    beta_group.slice(0) = Big_beta_response.rows(0, y.n_cols - 1).t();
    arma::mat beta_transform(p, (K + 1) * y.n_cols, arma::fill::zeros);
    beta_transform.cols(0, K) = arma::reshape(new_mat_group.slice(0).col(0), p, K + 1);
    arma::uvec input2 = arma::regspace<arma::uvec>(1, y.n_cols + 1 + K);
    arma::uvec multiple_of_K = modulo(input2, K + 1);
    arma::uvec II2 = input2.elem(arma::find(multiple_of_K == 0));
    int e2 = II2(0);
    // for (c_count2 in 2:dim(y)[2]) {
    //   beta_transform[,c((e2+1):(c_count2*(1+K)))]<-matrix(new.mat_group[,c_count2,1],ncol = (K+1), nrow = p)
    //   e2=II2[c_count2]
    // }


  //     norm_res<-((apply(beta_transform,c(1),twonorm)))
  //     coef.term1<- pmax(1-  (gg[1]) /rho/(norm_res),0)
  //     N_E1<-scale(t(beta_transform),center = FALSE,scale = 1/coef.term1)
  //     N_E1<-t(N_E1)
  //     beta_transform1<-matrix(0,p+p*K,dim(y)[2])
  //     beta_transform1[,1]<-as.vector(N_E1[,c(1:(K+1))])
  //     input3<-1:(dim(y)[2]*(1+K))
  //     multiple_of_K = (input3 %% (K+1)) == 0
  //     II3<-input3[multiple_of_K]
  //     e3=II3[-length(II3)][1]

  //     for (c_count3 in 2:dim(y)[2]) {
  //       beta_transform1[,c_count3]<-as.vector(N_E1[,c((e3+1):((K+1)*c_count3) )])
  //       e3=II3[c_count3]
  //     }

  //     N_E[[1]]<-(t(beta_transform1))
  //     e=II[-length(II)][1]

  //     for (c_count in 2:dim(C)[1]) {
  //       new.mat_group[,,c_count]<-t( (new.mat[c((e+1):( c_count*dim(y)[2]) ),]) )
  //       beta.group[,,c_count]<-t(Big_beta_respone[c((e+1):( c_count*dim(y)[2]) ),])
  //       beta_transform<-matrix(0,p,(K+1)*dim(y)[2])
  //       beta_transform[,c(1:(1+K))]<-matrix(new.mat_group[,1,c_count],ncol = (K+1), nrow = p)
  //       input2<-1:(dim(y)[2]*(1+K))
  //       multiple_of_K = (input2 %% (K+1)) == 0
  //       II2<-input2[multiple_of_K]
  //       e2=II2[-length(II2)][1]

  //       for (c_count2 in 2:dim(y)[2]) {
  //         beta_transform[,c((e2+1):(c_count2*(1+K)))]<-matrix(new.mat_group[,c_count2,c_count],ncol = (K+1), nrow = p)
  //         e2=II2[c_count2]
  //       }

  //       norm_res<-((apply(beta_transform,c(1),twonorm)))
  //       coef.term1<- pmax(1-  (gg[1]) /rho/(norm_res),0)
  //       N_E1<-scale(t(beta_transform),center = FALSE,scale = 1/coef.term1)
  //       N_E1<-t(N_E1)
  //       beta_transform1<-matrix(0,p+p*K,dim(y)[2])
  //       beta_transform1[,1]<-as.vector(N_E1[,c(1:(K+1))])
  //       input3<-1:(dim(y)[2]*(1+K))
  //       multiple_of_K = (input3 %% (K+1)) == 0
  //       II3<-input3[multiple_of_K]
  //       e3=II3[-length(II3)][1]

  //       for (c_count3 in 2:dim(y)[2]) {
  //         beta_transform1[,c_count3]<-as.vector(N_E1[,c((e3+1):((K+1)*c_count3) )])
  //         e3=II3[c_count3]
  //       }

  //       N_E[[c_count]]<-(t( (beta_transform1)) )
  //       e=II[c_count]
  //     } //end of for (c_count in 2:dim(C)[1])

  //     N_beta.group<-apply(beta.group, 3, twonorm)
  //     E[c(1:dim(C)[2]),]<-N_E[[1]]
  //     c_count<-2
  //     e=II[-length(II)][1]

  //     for (c_count in 2:dim(C)[1]) {
  //       E[c((e+1):( c_count*dim(y)[2]) ),]<-N_E[[c_count]]
  //       e=II[c_count]
  //     }

  //     H<-H+Big_beta_respone-E
  //     obj<-c(obj, obj)
  //     v.diff<-sum((-rho*(V-V_old))^2,na.rm = TRUE)
  //     q.diff<-sum((-rho*(Q-Q_old))^2,na.rm = TRUE)
  //     e.diff<-sum((-rho*(E-E_old))^2,na.rm = TRUE)
  //     ee.diff<-sum((-rho*(EE-EE_old))^2,na.rm = TRUE)
  //     s <- sqrt(v.diff+q.diff+e.diff+ee.diff)
  //     v.diff1<-sum(v.diff1)
  //     q.diff1<-sum(q.diff1)
  //     e.diff1<-sum(((Big_beta_respone-E))^2,na.rm = TRUE)
  //     ee.diff1<-sum(ee.diff1)#sum(((beta_hat-EE))^2,na.rm = TRUE)
  //     r <- sqrt(v.diff1+q.diff1+e.diff1+ee.diff1)
  //     res_dual<-s
  //     res_pri<-r
  //     e.primal <- sqrt(length(Big_beta11)+2*length(beta_hat) + length(Big_beta_respone) ) * e.abs + e.rel * max(twonorm(c((Big_beta11),(beta_hat),(beta_hat),(Big_beta_respone) )), twonorm(-c((V),(Q),(E),(EE) )))
  //     e.dual <-  sqrt(length(Big_beta11)+2*length(beta_hat)+length(Big_beta_respone) ) * e.abs + e.rel * twonorm((c((O),(P),(H),(HH) )))
  //     V_old <- V
  //     Q_old <- Q
  //     E_old<-E
  //     EE_old<-EE

  //     if( res_pri > 10*res_dual ){
  //       rho<- 2*rho
  //     }else if(res_pri*10 < res_dual ){
  //       rho<- rho/2
  //     }

  //     if(my_print==T){
  //       print(c(res_dual,e.dual,res_pri,e.primal))}

  //     if (res_pri <= e.primal && res_dual <= e.dual){
  //       # Update convergence message
  //       print(c("Convergence reached after  iterations",(i)))
  //       converge=T
  //       break
  //     }
    } // iteration; end of for (i in 2:max_it)



  // res_val<-t(I)%*%(E)
  // for (jj in 1:dim(y)[2]) {
  //   group<-(t(G)%*%t((V[,,jj]) )   )
  //   group1<-group[1,]; group2<-t(group[-1,])
  //   new_group=matrix(0,p,(K+1))
  //   new_group[,1]<-group1; new_group[,-1]<-group2
  //   new_g_theta<-as.vector(new_group)
  //   finB1<- as.vector(beta_hat[c(1:p),jj])*(new_g_theta[c(1:p)]!=0)*(as.vector((Q[,1,jj] ))!=0)
  //   finB2<- as.vector(beta_hat[-c(1:p),jj])*(new_g_theta[-c(1:p)]!=0)*(as.vector((Q[,-1,jj] ))!=0)
  //   beta_hat1<- matrix(c(finB1,finB2), ncol = (K+1), nrow = p )
  //   beta[,jj]<- beta_hat1[,1]#main_beta[,1]#
  //   theta[,,jj]<-(beta_hat1[,-1])
  //   beta_hat[,jj]<-c(c(beta_hat1[,1],as.vector(theta[,,jj])))
  // }

  arma::mat y_hat = Rcpp::as<arma::mat>(model_p(beta0.t(), theta0, beta_hat, theta, W_hat, Z));

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta0") = beta0,
    Rcpp::Named("theta0") = theta0,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("converge") = converge,
    Rcpp::Named("obj") = obj,
    Rcpp::Named("V") = V,
    Rcpp::Named("Q") = Q,
    Rcpp::Named("O") = O,
    Rcpp::Named("P") = P,
    Rcpp::Named("E") = E,
    Rcpp::Named("H") = H,
    Rcpp::Named("EE") = EE,
    Rcpp::Named("HH") = HH,
    Rcpp::Named("beta_hat") = beta_hat,
    Rcpp::Named("y_hat") = y_hat
  );
  out.attr("class") = "admm.MADMMplasso";
  return out;
}
