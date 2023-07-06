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
  const arma::mat beta_hat,
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
  const int K,
  const double e_abs,
  const double e_rel,
  const double alpha,
  const arma::vec lambda,
  const double alph,
  const Rcpp::List svd_w,
  const Rcpp::List tree,
  const Rcpp::List invmat,
  const arma::cube V,
  const arma::cube Q,
  const arma::mat E,
  const arma::cube EE,
  const arma::cube O,
  const arma::cube P,
  const arma::mat H,
  const arma::cube HH,
  const arma::vec gg,
  const bool my_print = true
) {
  Rcpp::List TT = tree;
  arma::sp_mat C = TT["Tree"]; // FIXME: sp_mat is for > 100x100 matrices. Use "mat" if this is not the case
  arma::vec CW = TT["Tw"];
  arma::mat svd_w_tu = Rcpp::as<arma::mat>(svd_w["u"]);
  arma::mat svd_w_tv = Rcpp::as<arma::mat>(svd_w["v"]);
  int D = y.n_cols;

  // for response groups #######################################################
  arma::ivec input = Rcpp::seq_len(D * C.n_rows);
  arma::ivec multiple_of_D = multiples_of(input, D);

  arma::mat I = arma::zeros<arma::mat>(C.n_rows * D, D);
  // arma::ivec II = input.elem(multiple_of_D);  // FIXME: cast multiple_of_D as uvec?
  //   II<-input[multiple_of_D]
  //   diag(I[c(1:dim(y)[2] ),])<-C[1,]*(CW[1])
  //   c_count<-2
  //   for (e in II[-length(II)]) {
  //     diag(I[c((e+1):( c_count*dim(y)[2]) ),]) <-C[c_count,]*(CW[c_count])
  //     c_count= 1+c_count
  //   }
  //   new_I=diag(t(I)%*%I)

  //   ###### overlapping group fro covariates########

  //   G=matrix(0,2*(1+K),(1+K))
  //   diag_G<-matrix(0,(K+1),(K+1))
  //   diag(diag_G)<-1
  //   for (i in 1:(K+1)) {
  //     G[i,]<-diag_G[i,]
  //   }
  //   for (i in (K+3):(2*(K+1))) {
  //     G[i,]<-diag_G[(i-(K+1)),]
  //   }
  //   V_old<-V;Q_old<-Q;E_old<-E;EE_old<-EE
  //   res_pri=0;res_dual=0
  arma::vec obj;
  //   SVD_D<-Diagonal(x=svd.w$d)
  //   R_svd<-(svd.w$u%*%SVD_D)/N

  //   rho=rho1
  //   Big_beta11<-V
  //   for (i in 2:max_it) {
  //     r_current = (y-model_intercept(beta0,theta0,beta=beta_hat, theta, X=W_hat, Z))
  //     b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
  //     beta0<-b$beta0
  //     theta0<-b$theta0
  //     new_y<- y-( matrix(1,N)%*%beta0+Z%*%((theta0)))
  //     XtY <- crossprod((W_hat), (new_y))
  //     main_beta<-array(0,c(p,K+1,D))
  //     res_val<-rho*(t(I)%*%(E)-(t(I)%*%(H)))
  //     v.diff1<-matrix(0,D); q.diff1<-matrix(0,D); ee.diff1<-matrix(0,D)

  //     new_G<-matrix(0,(p+p*K))
  //     new_G[c(1:p)]<-1;new_G[-c(1:p)]<-2
  //     new_G[c(1:p)]<-rho*(1+new_G[c(1:p)]);new_G[-c(1:p)]<-rho*(1+new_G[-c(1:p)])

  //     invmat<-list() #denominator of the beta estimates
  //     for (rr in 1:D) {
  //       DD1<-rho*(new_I[rr]+1)
  //       DD2<-new_G+DD1
  //       invmat[[rr]] <-DD2# Matrix::chol2inv( Matrix::chol(new_sparse) )
  //     }

  //     for (jj in 1:D) {
  //       group<-(rho)*(t(G)%*%t(V[,,jj])-t(G)%*%t(O[,,jj])   )
  //       group1<-group[1,]; group2<-t(group[-1,])
  //       new_group=matrix(0,p,(K+1))
  //       new_group[,1]<-group1; new_group[,-1]<-group2
  //       my_beta_jj<-XtY[,jj]/N  +as.vector(new_group)+as.vector(res_val[jj,])+as.vector(rho*(Q[,,jj]-P[,,jj] ))+as.vector(rho*(EE[,,jj]-HH[,,jj] ))
  //       my_beta_jj<-matrix(my_beta_jj,ncol = 1)
  //       DD3=Diagonal(x=1/invmat[[jj]])
  //       part_z<-DD3%*%t(W_hat )
  //       part_y<-DD3%*%my_beta_jj
  //       beta_hat_j<- solve(solve(R_svd)+(svd.w$tv)%*%part_z)
  //       beta_hat_j<-beta_hat_j%*%((svd.w$tv)%*%part_y)
  //       beta_hat_j<-part_z%*%beta_hat_j
  //       beta_hat_JJ<-part_y-beta_hat_j
  //       beta_hat_JJ<-matrix(beta_hat_JJ,ncol = 1)
  //       beta_hat[,jj]<-beta_hat_JJ
  //       beta_hat1<-matrix(beta_hat_JJ,p,(1+K))#main_beta[,,jj]#matrix(0,p,(K+1))
  //       b_hat<-alph*beta_hat1+(1-alph)*Q[,,jj]
  //       Q[,1,jj]<-b_hat[,1]+(P[,1,jj])
  //       new.mat<- b_hat[,-1] +P[,-1,jj]
  //       Q[,-1,jj]<- sign(new.mat)*pmax(abs(new.mat)-((alpha*lambda[jj])/(rho)),0)
  //       b_hat<-alph*beta_hat1+(1-alph)*EE[,,jj]
  //       new.mat<- b_hat +HH[,,jj]
  //       row.norm1<- sqrt(apply(new.mat^2,1,sum,na.rm = T))
  //       coef.term1<- pmax(1-  (gg[2]) /rho/(row.norm1),0)
  //     ee1<-scale(t(as.matrix(new.mat)),center = FALSE,scale = 1/coef.term1)
  //     EE[,,jj]<-  t(ee1)

  //       Big_beta<-t(tcrossprod(G,(beta_hat1)) )
  //       Big_beta11[,,jj]<-Big_beta
  //       Big_beta1<-alph*Big_beta+(1-alph)*V[,,jj]

  //       #Now we have the main part.
  //       new.mat<- Big_beta1+O[,,jj]
  //       new.mat1<-new.mat[,c(1:(K+1))];new.mat2<-new.mat[,-c(1:(K+1))]
  //       row.norm1<- sqrt(apply(new.mat1^2,1,sum,na.rm = T))
  //       row.norm2<- sqrt(apply(new.mat2^2,1,sum,na.rm = T))

  //       coef.term1<- pmax(1-( (1-alpha)*lambda[jj] )/(rho)/(row.norm1),0)
  //       coef.term2<- pmax(1-( (1-alpha)*lambda[jj] )/(rho)/(row.norm2),0)
  //       N_V1<-scale(t( new.mat1 ),center = FALSE,scale = 1/coef.term1)
  //       N_V2<-scale(t( new.mat2),center = FALSE,scale = 1/coef.term2)

  //       V[,,jj]= cbind(t(N_V1),t(N_V2))

  //       P[,,jj]<- P[,,jj]+beta_hat1-Q[,,jj]
  //       HH[,,jj]<-HH[,,jj]+beta_hat1-EE[,,jj]
  //       O[,,jj]<-O[,,jj]+Big_beta-V[,,jj]

  //       v.diff1[jj]<-sum(((Big_beta-V[,,jj]))^2,na.rm = TRUE)
  //       q.diff1[jj]<-sum(((beta_hat1-Q[,,jj]))^2,na.rm = TRUE)
  //       ee.diff1[jj]<-sum(((beta_hat1-EE[,,jj]))^2,na.rm = TRUE)
  //     }

  //     Big_beta_respone<-((I)%*%t( beta_hat ))
  //     b_hat_response<-alph*Big_beta_respone+(1-alph)*E
  //     new.mat<- b_hat_response +H

  //     new.mat_group<-(array(NA,c(p+p*K,dim(y)[2],dim(C)[1])))
  //     beta.group<-(array(NA,c(p+p*K,dim(y)[2],dim(C)[1])))
  //     N_E<-list()
  //     new.mat_group[,,1]<-t( (new.mat[c(1:dim(y)[2]),] ))
  //     beta.group[,,1]<-t( (Big_beta_respone[c(1:dim(y)[2]),]))

  //     beta_transform<-matrix(0,p,(K+1)*dim(y)[2])
  //     beta_transform[,c(1:(1+K))]<-matrix(new.mat_group[,1,1],ncol = (K+1), nrow = p)
  //     input2<-1:(dim(y)[2]*(1+K))
  //     multiple_of_K = (input2 %% (K+1)) == 0
  //     II2<-input2[multiple_of_K]
  //     e2=II2[-length(II2)][1]

  //     for (c_count2 in 2:dim(y)[2]) {
  //       beta_transform[,c((e2+1):(c_count2*(1+K)))]<-matrix(new.mat_group[,c_count2,1],ncol = (K+1), nrow = p)
  //       e2=II2[c_count2]
  //     }

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
  //     }

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
    bool converge = false;

  //   }### iteration



  //   res_val<-t(I)%*%(E)
  //   for (jj in 1:dim(y)[2]) {
  //     group<-(t(G)%*%t((V[,,jj]) )   )

  //     group1<-group[1,]; group2<-t(group[-1,])
  //     new_group=matrix(0,p,(K+1))
  //     new_group[,1]<-group1; new_group[,-1]<-group2
  //     new_g_theta<-as.vector(new_group)

  //     finB1<- as.vector(beta_hat[c(1:p),jj])*(new_g_theta[c(1:p)]!=0)*(as.vector((Q[,1,jj] ))!=0)
  //     finB2<- as.vector(beta_hat[-c(1:p),jj])*(new_g_theta[-c(1:p)]!=0)*(as.vector((Q[,-1,jj] ))!=0)

  //     beta_hat1<- matrix(c(finB1,finB2), ncol = (K+1), nrow = p )
  //     beta[,jj]<- beta_hat1[,1]#main_beta[,1]#
  //     theta[,,jj]<-(beta_hat1[,-1])
  //     beta_hat[,jj]<-c(c(beta_hat1[,1],as.vector(theta[,,jj])))
  //   }
  // Rcpp::Environment MADMMplasso = Rcpp::Environment::namespace_env("MADMMplasso"); // TODO: necessary?
  // Rcpp::Function model_p = MADMMplasso["model_p"]; // TODO: necessary?
  // Rcpp::NumericMatrix y_hat = model_p(beta0, theta0, beta_hat, theta, W_hat, Z); // FIXME: output as arma::mat
  Rcpp::NumericMatrix y_hat; // TEMP: placeholder

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
