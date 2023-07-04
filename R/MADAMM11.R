#' Fit the ADMM part of  model for a given lambda vale
#' @param X  n by p matrix of predictors
#' @param Z n by nz matrix of modifying variables. The elements of z
#' may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical varables should be coded by 0-1 dummy variables: for a k-level
#' variable, one can use either k or k-1  dummy variables.
#' @param beta0 TODO: fill paramater description
#' @param theta0 TODO: fill paramater description
#' @param beta TODO: fill paramater description
#' @param beta_hat TODO: fill paramater description
#' @param theta TODO: fill paramater description
#' @param rho1 TODO: fill paramater description
#' @param max_it TODO: fill paramater description
#' @param W_hat TODO: fill paramater description
#' @param XtY TODO: fill paramater description
#' @param y TODO: fill paramater description
#' @param N TODO: fill paramater description
#' @param p TODO: fill paramater description
#' @param K TODO: fill paramater description
#' @param e.abs TODO: fill paramater description
#' @param e.rel TODO: fill paramater description
#' @param alpha TODO: fill paramater description
#' @param lambda TODO: fill paramater description
#' @param alph TODO: fill paramater description
#' @param svd.w TODO: fill paramater description
#' @param tree TODO: fill paramater description
#' @param my_print TODO: fill paramater description
#' @param invmat TODO: fill paramater description
#' @param V TODO: fill paramater description
#' @param Q TODO: fill paramater description
#' @param E TODO: fill paramater description
#' @param EE TODO: fill paramater description
#' @param O TODO: fill paramater description
#' @param P TODO: fill paramater description
#' @param H TODO: fill paramater description
#' @param HH TODO: fill paramater description
#' @param cv TODO: fill paramater description
#' @param gg TODO: fill paramater description
#' @return  predicted values for the ADMM part
#' @description TODO: add description



#' @export
admm.MADMMplasso<-function(beta0,theta0,beta,beta_hat,theta,rho1,X,Z,max_it,W_hat,XtY,y,N,p,K,e.abs, e.rel,alpha,lambda,alph,svd.w,tree,my_print=TRUE,invmat,V,Q,E,EE,O,P,H,HH,cv=cv,gg=0.2){

  TT<-tree

  C<-TT$Tree
  CW<-TT$Tw
  svd.w$tu<- t(svd.w$u)
  svd.w$tv<- t(svd.w$v)
  D=dim(y)[2]







  ### for response groups ###############################################################


  input<-1:(dim(y)[2]*nrow(C))
  multiple_of_D = (input %% dim(y)[2]) == 0




  I<-matrix(0,nrow = nrow(C)*dim(y)[2],ncol = dim(y)[2])
  II<-input[multiple_of_D]
  diag(I[c(1:dim(y)[2] ),])<-C[1,]*(CW[1])

  c_count<-2
  for (e in II[-length(II)]) {
    #	print(e+1)
    #print(diag(I[c((e+1):( c_count*dim(y)[2]) ),]))
    diag(I[c((e+1):( c_count*dim(y)[2]) ),]) <-C[c_count,]*(CW[c_count])
    c_count= 1+c_count
  }


  new_I=diag(t(I)%*%I)


  #########################################################################
  ###### overlapping group fro covariates########


  G=matrix(0,2*(1+K),(1+K))
  diag_G<-matrix(0,(K+1),(K+1))
  diag(diag_G)<-1
  for (i in 1:(K+1)) {
    G[i,]<-diag_G[i,]
  }
  for (i in (K+3):(2*(K+1))) {
    G[i,]<-diag_G[(i-(K+1)),]
  }

  ####################################################################

  # auxiliary variables for the L1 norm####



  #new_I=diag(t(I)%*%I)
  ################################################

  V_old<-V;Q_old<-Q;E_old<-E;EE_old<-EE
  res_pri=0;res_dual=0
  obj<-c()





  SVD_D<-Diagonal(x=svd.w$d)
  R_svd<-(svd.w$u%*%SVD_D)/N


  rho=rho1
  #
  #XtY <- crossprod(my_W_hat, y)
  Big_beta11<-V
  for (i in 2:max_it) {






    r_current = (y-model_intercept(beta0,theta0,beta=beta_hat, theta, X=W_hat, Z))
    b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
    beta0<-b$beta0
    theta0<-b$theta0

    new_y<- y-( matrix(1,N)%*%beta0+Z%*%((theta0)))

    XtY <- crossprod((W_hat), (new_y))


    #W_hat.s<-scale(W_hat)
    #r_current = y-model_intercept(beta=beta_hat, theta, X=W_hat, Z)


    #factor(my_W_hat,rho)
    main_beta<-array(0,c(p,K+1,D))
    # r_current = y-model(beta0, theta0, beta, theta, X, Z)

    # print(theta1)
    #	group=matrix(0,p,(K+1))
    res_val<-rho*(t(I)%*%(E)-(t(I)%*%(H)))
    #res_val1<-rho*(EE-HH)




    v.diff1<-matrix(0,D); q.diff1<-matrix(0,D); ee.diff1<-matrix(0,D)

    new_G<-matrix(0,(p+p*K))
    new_G[c(1:p)]<-1;new_G[-c(1:p)]<-2
    new_G[c(1:p)]<-rho*(1+new_G[c(1:p)]);new_G[-c(1:p)]<-rho*(1+new_G[-c(1:p)])

    invmat<-list() #denominator of the beta estimates
    for (rr in 1:D) {

      DD1<-rho*(new_I[rr]+1)

      #DD1<-rho1*(obv_beta_matrix)

      DD2<-new_G+DD1

      #DD2[c(1:p)]<-DD2[c(1:p)]+DD1






      #print(rbind(part4[3,]-as.matrix(part4[3,]) ))
      invmat[[rr]] <-DD2# Matrix::chol2inv( Matrix::chol(new_sparse) )
      #print(dim(invmat))
      # # #Matrix::chol2inv()


    }




    for (jj in 1:D) {



      group<-(rho)*(t(G)%*%t(V[,,jj])-t(G)%*%t(O[,,jj])   )
      #group1<-group[1,]+as.vector(res_val[jj,]); group2<-t(group[-1,])
      group1<-group[1,]; group2<-t(group[-1,])
      new_group=matrix(0,p,(K+1))

      new_group[,1]<-group1; new_group[,-1]<-group2


      #my_beta_jj<-XtY[,jj]/N  +as.vector(new_group) +as.vector(rho*(Q[,,jj]-P[,,jj] ))
      # my_beta_jj<-XtY[,jj]/N  +as.vector(new_group)+as.vector(res_val[jj,])+as.vector(res_val1[,jj])+as.vector(rho*(Q[,,jj]-P[,,jj] ))
      my_beta_jj<-XtY[,jj]/N  +as.vector(new_group)+as.vector(res_val[jj,])+as.vector(rho*(Q[,,jj]-P[,,jj] ))+as.vector(rho*(EE[,,jj]-HH[,,jj] ))

      my_beta_jj<-matrix(my_beta_jj,ncol = 1)


      DD3=Diagonal(x=1/invmat[[jj]])


      part_z<-DD3%*%t(W_hat )
      part_y<-DD3%*%my_beta_jj

      beta_hat_j<- solve(solve(R_svd)+(svd.w$tv)%*%part_z)
      beta_hat_j<-beta_hat_j%*%((svd.w$tv)%*%part_y)
      beta_hat_j<-part_z%*%beta_hat_j
      beta_hat_JJ<-part_y-beta_hat_j
      beta_hat_JJ<-matrix(beta_hat_JJ,ncol = 1)


      #beta_hat_JJ<-as.matrix(tcrossprod(invmat[[jj]],t((my_beta_jj))))

      # beta_hat[,jj]<-beta_hat_JJ
     # main_beta[,,jj]<-matrix(beta_hat_JJ,p,(1+K))
      beta_hat[,jj]<-beta_hat_JJ
      #main_beta[,-1,jj]<-matrix(beta_hat[,jj][-(1:p)],ncol = (K),nrow = p, byrow = T)
     # beta[,jj]<- main_beta[,1,jj]#beta_hat[c(1:p)]
      #beta_theta<-main_beta[,-1]
      #theta[,,jj]<-main_beta[,-1,jj]


      beta_hat1<-matrix(beta_hat_JJ,p,(1+K))#main_beta[,,jj]#matrix(0,p,(K+1))

      #beta_hat1[,1]<-beta; beta_hat1[,c(2:(K+1))]<-beta_theta
      #for (j in 1:p) {
      b_hat<-alph*beta_hat1+(1-alph)*Q[,,jj]
      Q[,1,jj]<-b_hat[,1]+(P[,1,jj])
      new.mat<- b_hat[,-1] +P[,-1,jj]
      Q[,-1,jj]<- sign(new.mat)*pmax(abs(new.mat)-((alpha*lambda[jj])/(rho)),0)
      # coef.term <- pmax(1-(lambda[1]/rho)/(row.norm) , 0)
      b_hat<-alph*beta_hat1+(1-alph)*EE[,,jj]
      new.mat<- b_hat +HH[,,jj]
      row.norm1<- sqrt(apply(new.mat^2,1,sum,na.rm = T))
      coef.term1<- pmax(1-  (gg[2]) /rho/(row.norm1),0)
    ee1<-scale(t(as.matrix(new.mat)),center = FALSE,scale = 1/coef.term1)
    #print(dim(ee1))
   # print(dim(EE[,,jj]))
    EE[,,jj]<-  t(ee1)
      #EE[,,jj]<- sign(new.mat)*pmax(abs(new.mat)-(( gg[1] )/(rho)),0)
      #row.norm1<- sqrt(apply(new.mat^2,1,sum,na.rm = T))
      # coef.term1<- pmax(1-( gg[1] )/(rho)/(row.norm1),0)
      # N_V1<-scale(t(new.mat),center = FALSE,scale = 1/coef.term1)
      # EE[,,jj]<-t(N_V1)


      Big_beta<-t(tcrossprod(G,(beta_hat1)) )
      Big_beta11[,,jj]<-Big_beta
      Big_beta1<-alph*Big_beta+(1-alph)*V[,,jj]


      #Now we have the main part.
      new.mat<- Big_beta1+O[,,jj]
      new.mat1<-new.mat[,c(1:(K+1))];new.mat2<-new.mat[,-c(1:(K+1))]
      row.norm1<- sqrt(apply(new.mat1^2,1,sum,na.rm = T))
      row.norm2<- sqrt(apply(new.mat2^2,1,sum,na.rm = T))


      coef.term1<- pmax(1-( (1-alpha)*lambda[jj] )/(rho)/(row.norm1),0)
      coef.term2<- pmax(1-( (1-alpha)*lambda[jj] )/(rho)/(row.norm2),0)
      #new.mat3<-matrix(0,p,(K+1));new.mat4<-matrix(0,p,(K+1))
      N_V1<-scale(t( new.mat1 ),center = FALSE,scale = 1/coef.term1)
      N_V2<-scale(t( new.mat2),center = FALSE,scale = 1/coef.term2)
      # N_V<-   lapply(seq_len(p),
      #           function(j) (  c(matrix(scale(new.mat1[j,], center = FALSE, scale = 1/coef.term1[j])  ),matrix(scale(new.mat2[j,], center = FALSE, scale = 1/coef.term2[j]) ) )  ) )



      V[,,jj]= cbind(t(N_V1),t(N_V2))


      #print(V)
      P[,,jj]<- P[,,jj]+beta_hat1-Q[,,jj]
      HH[,,jj]<-HH[,,jj]+beta_hat1-EE[,,jj]
      O[,,jj]<-O[,,jj]+Big_beta-V[,,jj]


      v.diff1[jj]<-sum(((Big_beta-V[,,jj]))^2,na.rm = TRUE)
      q.diff1[jj]<-sum(((beta_hat1-Q[,,jj]))^2,na.rm = TRUE)
      ee.diff1[jj]<-sum(((beta_hat1-EE[,,jj]))^2,na.rm = TRUE)


    }



    ############to estimate E  ##################
    ############


    # b_hat<-alph*beta_hat+(1-alph)*EE
    #Q[,1,jj]<-b_hat[,1]+(P[,1,jj])
    # new.mat<- b_hat +HH
    #EE<- sign(new.mat)*pmax(abs(new.mat)-(( gg[1] )/(rho)),0)

    #b_hat<-alph*beta_hat+(1-alph)*EE
    #Q[,1,jj]<-b_hat[,1]+(P[,1,jj])
    #new.mat<- b_hat +HH
    # EE<- sign(new.mat)*pmax(abs(new.mat)-( (gg[2])/rho),0)


    # row.norm1<- (sqrt(apply((new.mat)^2,1,sum,na.rm = T)))



    # coef.term1<- pmax(1-((1-alpha)*min(lambda) )/(rho)/(row.norm1),0)

    #new.mat3<-matrix(0,p,(K+1));new.mat4<-matrix(0,p,(K+1))
    # EE<-t(scale(t(new.mat),center = FALSE,scale = 1/coef.term1))




    Big_beta_respone<-((I)%*%t( beta_hat ))
    b_hat_response<-alph*Big_beta_respone+(1-alph)*E
    #Q[,1]<-b_hat[,1]-(P[,1])/rho
    new.mat<- b_hat_response +H

    #row.norm1<- (sqrt(apply(new.mat^2,1,sum,na.rm = T)))
    #coef.term1<- pmax(1- (((1-alpha)* min(lambda) )/rho)/(row.norm1),0)
    #E<-t(scale(t(new.mat),center = FALSE,scale = 1/coef.term1))
    #print(dim(E))
    #print(dim(H))
    #                  r
    #
    new.mat_group<-(array(NA,c(p+p*K,dim(y)[2],dim(C)[1])))
    beta.group<-(array(NA,c(p+p*K,dim(y)[2],dim(C)[1])))
    N_E<-list()
    #I<-matrix(0,nrow = nrow(C)*dim(y)[2],ncol = dim(y)[2])
    II<-input[multiple_of_D]
    new.mat_group[,,1]<-t( (new.mat[c(1:dim(y)[2]),] ))
    beta.group[,,1]<-t( (Big_beta_respone[c(1:dim(y)[2]),]))


    beta_transform<-matrix(0,p,(K+1)*dim(y)[2])
    beta_transform[,c(1:(1+K))]<-matrix(new.mat_group[,1,1],ncol = (K+1), nrow = p)
    input2<-1:(dim(y)[2]*(1+K))
    multiple_of_K = (input2 %% (K+1)) == 0
    II2<-input2[multiple_of_K]
    e2=II2[-length(II2)][1]

    for (c_count2 in 2:dim(y)[2]) {

      beta_transform[,c((e2+1):(c_count2*(1+K)))]<-matrix(new.mat_group[,c_count2,1],ncol = (K+1), nrow = p)



      e2=II2[c_count2]

    }


    norm_res<-((apply(beta_transform,c(1),twonorm)))
    coef.term1<- pmax(1-  (gg[1]) /rho/(norm_res),0)

    N_E1<-scale(t(beta_transform),center = FALSE,scale = 1/coef.term1)

    N_E1<-t(N_E1)
    beta_transform1<-matrix(0,p+p*K,dim(y)[2])
    beta_transform1[,1]<-as.vector(N_E1[,c(1:(K+1))])

    input3<-1:(dim(y)[2]*(1+K))
    multiple_of_K = (input3 %% (K+1)) == 0
    II3<-input3[multiple_of_K]
    e3=II3[-length(II3)][1]

    for (c_count3 in 2:dim(y)[2]) {

      beta_transform1[,c_count3]<-as.vector(N_E1[,c((e3+1):((K+1)*c_count3) )])



      e3=II3[c_count3]

    }

    #print(beta_transform1)
    N_E[[1]]<-(t(beta_transform1))


    e=II[-length(II)][1]
    for (c_count in 2:dim(C)[1]) {

      #for (e in II[-length(II)]) {
      new.mat_group[,,c_count]<-t( (new.mat[c((e+1):( c_count*dim(y)[2]) ),]) )
      beta.group[,,c_count]<-t(Big_beta_respone[c((e+1):( c_count*dim(y)[2]) ),])
      # c_count= 1+c_count
      # }
      #e=II[-length(II)][1+1]
      #e=II[c_count]




      beta_transform<-matrix(0,p,(K+1)*dim(y)[2])
      beta_transform[,c(1:(1+K))]<-matrix(new.mat_group[,1,c_count],ncol = (K+1), nrow = p)
      input2<-1:(dim(y)[2]*(1+K))
      multiple_of_K = (input2 %% (K+1)) == 0
      II2<-input2[multiple_of_K]
      e2=II2[-length(II2)][1]

      for (c_count2 in 2:dim(y)[2]) {

        beta_transform[,c((e2+1):(c_count2*(1+K)))]<-matrix(new.mat_group[,c_count2,c_count],ncol = (K+1), nrow = p)



        e2=II2[c_count2]

      }


      norm_res<-((apply(beta_transform,c(1),twonorm)))
      coef.term1<- pmax(1-  (gg[1]) /rho/(norm_res),0)

      N_E1<-scale(t(beta_transform),center = FALSE,scale = 1/coef.term1)

      N_E1<-t(N_E1)
      beta_transform1<-matrix(0,p+p*K,dim(y)[2])
      beta_transform1[,1]<-as.vector(N_E1[,c(1:(K+1))])

      input3<-1:(dim(y)[2]*(1+K))
      multiple_of_K = (input3 %% (K+1)) == 0
      II3<-input3[multiple_of_K]
      e3=II3[-length(II3)][1]

      for (c_count3 in 2:dim(y)[2]) {

        beta_transform1[,c_count3]<-as.vector(N_E1[,c((e3+1):((K+1)*c_count3) )])



        e3=II3[c_count3]

      }

      #print(beta_transform1)
      N_E[[c_count]]<-(t( (beta_transform1)) )



      e=II[c_count]

    }

    # new.mat_group<-array(NA,c(p+p*K,dim(y)[2],dim(C)[1]))
    # beta.group<-array(NA,c(p+p*K,dim(y)[2],dim(C)[1]))
    # #I<-matrix(0,nrow = nrow(C)*dim(y)[2],ncol = dim(y)[2])
    # II<-input[multiple_of_D]
    # new.mat_group[,,1]<-t(new.mat[c(1:dim(y)[2]),])
    # beta.group[,,1]<-t(Big_beta_respone[c(1:dim(y)[2]),])
    # c_count<-2
    # e=II[-length(II)][1]
    # for (c_count in 2:dim(C)[1]) {
    #
    #   #for (e in II[-length(II)]) {
    #   new.mat_group[,,c_count]<-t(new.mat[c((e+1):( c_count*dim(y)[2]) ),])
    #   beta.group[,,c_count]<-t(Big_beta_respone[c((e+1):( c_count*dim(y)[2]) ),])
    #   # c_count= 1+c_count
    #   # }
    #   #e=II[-length(II)][1+1]
    #   e=II[c_count]
    #
    # }
    #for (cc in 1:dim(C)[1]) {
     # N_E<-   lapply(seq_len(dim(C)[1]),
     #                function(g){
     #                  row.norm1<- (sqrt(apply(new.mat_group[,,g]^2,1,sum,na.rm = T)))
     #
     #                  coef.term1<- pmax(1-  (gg[1]) /rho/(row.norm1),0)
     #                  N_E1<-scale(t(new.mat_group[,,g]),center = FALSE,scale = 1/coef.term1)
     #                 return(N_E1)        })
    #norm_res<-((apply(new.mat_group,c(1),twonorm)))
    #coef.term1<- pmax(1-  (gg[1]) /rho/(norm_res),0)
    #print(dim(coef.term1))
    #coef.term1<-matrix(1,p)%*%coef.term1


   # N_E<-   lapply(seq_len(dim(C)[1]),
    #               function(g){




    #                 norm_res<-((apply(new.mat_group[,,g],c(1),twonorm)))
    #                 coef.term1<- pmax(1-  (gg[1]) /rho/(norm_res),0)

    #                 N_E1<-scale(t(new.mat_group[,,g]),center = FALSE,scale = 1/coef.term1)

                     #print(beta_transform1)

                     #print(N_E1)
 #                    return(N_E1)        })



    N_beta.group<-apply(beta.group, 3, twonorm)

    E[c(1:dim(C)[2]),]<-N_E[[1]]


    c_count<-2
    e=II[-length(II)][1]
    #ee=1
    for (c_count in 2:dim(C)[1]) {

      #for (e in II[-length(II)]) {
      E[c((e+1):( c_count*dim(y)[2]) ),]<-N_E[[c_count]]
      #c_count= 1+c_count
      # }
      #ee=ee+1
      e=II[c_count]

    }

    # coef.term1<- pmax(1-((1-alpha)*lambda)/(rho)/(row.norm1),0)

    #new.mat3<-matrix(0,p,(K+1));new.mat4<-matrix(0,p,(K+1))
    # N_E<-scale(t(new.mat1),center = FALSE,scale = 1/coef.term1)


    #}



    #E<- sign(new.mat)*pmax(abs(new.mat)-((lambda)/(rho)),0)

    H<-H+Big_beta_respone-E
    ################################################


    # print(V)
    #print(Big_beta_respone)
    #obj1<-objective(beta0,theta0,
    #               beta,theta,X,Z,
    #            y,alpha,lambda=lambda,p,N,IB=sum(unlist(N_beta.group)),W_hat,beta1=beta_hat11)
    #print(obj1)
    obj<-c(obj, obj)

    # Calculate residuals for iteration t
    # print(Big_beta)
    # print(V)
    # r <- cbind(Big_beta,beta_hat1)-cbind(V,Q)

    v.diff<-sum((-rho*(V-V_old))^2,na.rm = TRUE)

    q.diff<-sum((-rho*(Q-Q_old))^2,na.rm = TRUE)
    e.diff<-sum((-rho*(E-E_old))^2,na.rm = TRUE)
    ee.diff<-sum((-rho*(EE-EE_old))^2,na.rm = TRUE)

    s <- sqrt(v.diff+q.diff+e.diff+ee.diff)

    v.diff1<-sum(v.diff1)
    q.diff1<-sum(q.diff1)
    e.diff1<-sum(((Big_beta_respone-E))^2,na.rm = TRUE)
    ee.diff1<-sum(ee.diff1)#sum(((beta_hat-EE))^2,na.rm = TRUE)
    r <- sqrt(v.diff1+q.diff1+e.diff1+ee.diff1)

    #print(i)
    #print(V)
    #print(V_old)
    # print(O)

    res_dual<-s
    res_pri<-r
    #print(c(res_dual,res_pri))

    e.primal <- sqrt(length(Big_beta11)+2*length(beta_hat) + length(Big_beta_respone) ) * e.abs + e.rel * max(twonorm(c((Big_beta11),(beta_hat),(beta_hat),(Big_beta_respone) )), twonorm(-c((V),(Q),(E),(EE) )))

    #print(c(twonorm(cbind(Big_beta)), twonorm(-cbind(V))))
    e.dual <-  sqrt(length(Big_beta11)+2*length(beta_hat)+length(Big_beta_respone) ) * e.abs + e.rel * twonorm((c((O),(P),(H),(HH) )))
    V_old <- V
    Q_old <- Q
    E_old<-E
    EE_old<-EE

    if( res_pri > 10*res_dual ){
      rho<- 2*rho
    }else if(res_pri*10 < res_dual ){
      rho<- rho/2
    }

    if(my_print==T){
      print(c(res_dual,e.dual,res_pri,e.primal))}
    #print(c(res_dual,res_pri,sqrt(N*2*(1+K))*e.rel))
    if (res_pri <= e.primal && res_dual <= e.dual){
      # Remove excess beta and nll


      # Update convergence message
      print(c("Convergence reached after  iterations",(i)))
      converge=T
      break
    }
    converge = F

  }### iteration



  res_val<-t(I)%*%(E)
  #res_val1<-EE
  for (jj in 1:dim(y)[2]) {
    group<-(t(G)%*%t((V[,,jj]) )   )


    # group1<-group[1,]+as.vector(res_val[jj,])+as.vector(res_val1[,jj]); group2<-t(group[-1,])
    group1<-group[1,]; group2<-t(group[-1,])
    #group1<-group[1,]; group2<-t(group[-1,])
    new_group=matrix(0,p,(K+1))
    new_group[,1]<-group1; new_group[,-1]<-group2
    new_g_theta<-as.vector(new_group)
    #group<-as.vector(new_group) +as.vector(Q)

    #finB1<- as.vector(beta_hat11[,jj])*(new_g_theta!=0)*(as.vector((Q[,,jj] ))!=0)


    # beta_hat111<- matrix((finB1), ncol = (K+1), nrow = p )


    #finB1<- as.vector(beta_hat[c(1:p),jj])*(new_g_theta[c(1:p)]!=0)*(as.vector((Q[,1,jj] ))!=0)*(as.vector((res_val[jj,c(1:p)] ))!=0)
    # finB1<- as.vector(beta_hat[c(1:p),jj])*(new_g_theta[c(1:p)]!=0)*(as.vector((Q[,1,jj] ))!=0)*(as.vector((EE[,jj] ))!=0)
    finB1<- as.vector(beta_hat[c(1:p),jj])*(new_g_theta[c(1:p)]!=0)*(as.vector((Q[,1,jj] ))!=0)
    #finB2<- as.vector(beta_hat[-c(1:p),jj])*(new_g_theta[-c(1:p)]!=0)*(as.vector((Q[,-1,jj] ))!=0)*(as.vector((res_val[jj,-c(1:p)] ))!=0)
    finB2<- as.vector(beta_hat[-c(1:p),jj])*(new_g_theta[-c(1:p)]!=0)*(as.vector((Q[,-1,jj] ))!=0)

    beta_hat1<- matrix(c(finB1,finB2), ncol = (K+1), nrow = p )
    #  # pinv(t(my_w)%*%my_w)
    #main_beta<-matrix(beta_hat,p,(1+K))
    beta[,jj]<- beta_hat1[,1]#main_beta[,1]#
    # print(c(beta_hat11[c(1:p),jj]))
    theta[,,jj]<-(beta_hat1[,-1])
    beta_hat[,jj]<-c(c(beta_hat1[,1],as.vector(theta[,,jj])))
  }
  # print(dim(W_hat)); print(dim(beta_hat11))
  y_hat<-model_p(beta0, theta0, beta=beta_hat, theta, X=W_hat, Z)

  out=list(beta0=beta0,theta0=theta0,beta=beta,theta=theta,converge=converge,obj=obj,V=V,Q=Q,O=O,P=P,E=E,H=H,EE=EE,HH=HH,beta_hat=beta_hat,y_hat=y_hat)
  class(out)="admm.MADMMplasso"

  return(out)


}


#' @title Fit a multi-response pliable lasso model over a path of regularization values
#' @description TODO: add description (This function fits a multi-response pliable lasso model over a path of regularization values?)
#' @param X  N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical varables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param y N by D matrix  of responses. The X and Z variables are centered in the function. We recommmend that X and Z also be standardized before the call
#' @param nlambda number of lambda values desired (default 50).
#' @param alpha mixing parameter- default 0.5
#' @param max_it maximum number of iterations in the ADMM algorithm for one lambda. Default 50000
#' @param maxgrid similar to nlambda
#' @param rho the Lagrange variable for the ADMM
#' @param e.abs absolute error for the admm
#' @param e.rel relative error for the admm
#' @param gg penalty term for the tree structure
#' @param my_lambda TODO: fill in
#' @param lambda_min TODO: fill in
#' @param max_it TODO: fill in
#' @param my_print TODO: fill in
#' @param alph TODO: fill in
#' @param tree TODO: fill in
#' @param cv TODO: fill in
#' @param parallel TODO: fill in
#' @param pal TODO: fill in
#' @param tol TODO: fill in
#' @param cl TODO: fill in
#' @return  predicted values for the MADMMplasso fit
#' @examples
#' # Train the model
#' # generate some data
#' set.seed(1235)
#' N = 100 ; p =50;nz=4; K=nz
#' X <- matrix(rnorm(n = N * p), nrow = N, ncol = p)
#' mx=colMeans(X)
#' sx=sqrt(apply(X,2,var))
#' X=scale(X,mx,sx)
#' X=matrix(as.numeric(X),N,p)
#' Z =matrix(rnorm(N*nz),N,nz)
#' mz=colMeans(Z)
#' sz=sqrt(apply(Z,2,var))
#' Z=scale(Z,mz,sz)
#' beta_1 <- rep(x = 0, times = p)
#' beta_2<-rep(x = 0, times = p)
#' beta_3<-rep(x = 0, times = p)
#' beta_4<-rep(x = 0, times = p)
#' beta_5<-rep(x = 0, times = p)
#' beta_6<-rep(x = 0, times = p)
#'
#'
#' beta_1[1:5] <- c(2, 2, 2, 2,2)
#' beta_2[1:5]<-c(2, 2, 2, 2,2)
#' beta_3[6:10]<-c(2, 2, 2, -2,-2)
#' beta_4[6:10] <- c(2, 2, 2, -2,-2)
#' beta_5[11:15] <- c(-2,  -2,-2, -2,-2)
#' beta_6[11:15] <- c(-2, -2, -2, -2,-2)
#'
#' Beta<-cbind(beta_1,beta_2,beta_3,beta_4,beta_5,beta_6)
#' colnames(Beta)<-c(1:6)
#'
#'
#' theta<-array(0,c(p,K,6))
#' theta[1,1,1]<-2;theta[3,2,1]<-2;theta[4,3,1]<- -2;theta[5,4,1]<- -2;
#' theta[1,1,2]<-2;theta[3,2,2]<-2;theta[4,3,2]<- -2;theta[5,4,2]<- -2;
#' theta[6,1,3]<-2;theta[8,2,3]<-2;theta[9,3,3]<- -2;theta[10,4,3]<- -2;
#' theta[6,1,4]<-2;theta[8,2,4]<-2;theta[9,3,4]<- -2;theta[10,4,4]<- -2;
#' theta[11,1,5]<-2;theta[13,2,5]<-2;theta[14,3,5]<- -2;theta[15,4,5]<- -2;
#' theta[11,1,6]<-2;theta[13,2,6]<-2;theta[14,3,6]<- -2;theta[15,4,6]<- -2
#'
#' library(MASS)
#' pliable = matrix(0,N,6)
#' for (e in 1:6) {
#' pliable[,e]<-	compute_pliable(X, Z, theta[,,e])
#' }
#'
#' esd<-diag(6)
#' e<-MASS::mvrnorm(N,mu=rep(0,6),Sigma=esd)
#' y_train<-X%*%Beta+pliable+e
#' y=y_train
#' #colnames(y)<-c(1:6)
#' colnames(y)<- c( paste("y",1:(ncol(y)),sep = "") )
#' TT=tree.parms(y)
#' plot(TT$h_clust)
#' gg1=matrix(0,2,2)
#' gg1[1,]<-c(0.02,0.02)
#' gg1[2,]<-c(0.02,0.02)
#'
#' nlambda = 1
#' e.abs=1E-4
#' e.rel=1E-2
#' alpha=.2
#' tol=1E-3
#' fit <- MADMMplasso(
#'   X, Z, y, alpha=alpha, my_lambda=matrix(rep(0.2,dim(y)[2]),1),
#'   lambda_min=0.001, max_it=5000, e.abs=e.abs, e.rel=e.rel, maxgrid=nlambda,
#'   nlambda=nlambda, rho=5, tree=TT, my_print=FALSE, alph=1, parallel=FALSE,
#'   pal=1, gg=gg1, tol=tol, cl=6
#' )
#' plot(fit)
#' @export
MADMMplasso<-function(X,Z,y,alpha,my_lambda=NULL,lambda_min=.001,max_it=50000,e.abs=1E-3,e.rel=1E-3,maxgrid,nlambda, rho=5,my_print=FALSE,alph=1.8,tree,cv=FALSE,parallel=TRUE,pal=0,gg=NULL,tol=1E-4,cl=4){

  N=nrow(X)
  #print(c(N,length(y)))

  p=ncol(X)
  K=ncol(Z)
  D=dim(y)[2]



  BETA0<-lapply(seq_len(nlambda),
                function(j)(matrix(0,nrow = (D))))

  BETA<-lapply(seq_len(nlambda),
               function(j)(as(matrix(0,nrow=p,ncol=(D) ),"sparseMatrix") ) )
  BETA_hat<-lapply(seq_len(nlambda),
                   function(j)( as(matrix(0,nrow=p+p*K,ncol=(D) ),"sparseMatrix") ) )


  THETA0<-lapply(seq_len(nlambda),
                 function(j)(matrix(0,nrow=K,ncol=(D)  )))

  THETA<-lapply(seq_len(nlambda),
                function(j)( as.sparse3Darray(array(0,c(p,K,(D)  ))) ))

  Y_HAT<-lapply(seq_len(nlambda),
                function(j)( matrix(0,nrow=N,ncol=(D) )  )  )
  #W<-compute_w(X, Z)
  #W_hat<-compute_w_hat(X, Z)
  # if(cv==T){
  # mx=colMeans(X)
  #
  # sx=sqrt(apply(X,2,var))
  # X=scale(X,mx,sx)
  # X=matrix(as.numeric(X),N,p)
  #
  # mz=colMeans(Z)
  # sz=sqrt(apply(Z,2,var))
  #
  # Z=scale(Z,mz,sz)
  # }
  #
  rat=lambda_min


  if(is.null(my_lambda)){
    lamda_new<-matrix(0,dim(y)[2])
    r = y

    lammax<-	lapply(seq_len(dim(y)[2]),
                    function(g){l_max<- max(abs(t(X)%*%(r- colMeans(r) ) )/length(r[,1]))/((1-alpha))
                    # l_max<-l_max/((1-alpha)*l_max/gg)
                    return( l_max)
                    })


    #lammax=max(abs(t(X)%*%r[,1])/length(r[,1]))/((1-alpha))
    # lam<-mean(unlist(lammax))
    #
    # lammax<-	lapply(seq_len(dim(y)[2]),
    #                 function(g){l_max<- lam
    #l_max<-l_max/((1-alpha)*l_max/gg)
    #                return( l_max)
    #              })
    #
    big_lambda<-lammax
    #print(big_lambda)

    lambda_i<-	lapply(seq_len(dim(y)[2]),
                      function(g){
                        lam_i<- exp(seq(log(big_lambda[[g]]),log(big_lambda[[g]]*rat),length=maxgrid))

                        return(lam_i)
                      }      )
   # gg1<-	lapply(seq_len(dim(y)[2]),
  #                  function(g){l_max<- max(abs(t(X)%*%(r- colMeans(r) ) )/length(r[,1]))
                    # l_max<-l_max/((1-alpha)*l_max/gg)
   #                return( l_max)
   #                 })
   # gg1= exp(seq(log(gg1[[1]]),log(gg1[[1]]*rat1),length=maxgrid))
    gg1= seq((gg[1,1]),(gg[1,2]),length=maxgrid)
    gg2= seq((gg[2,1]),(gg[2,2]),length=maxgrid)
    gg3=matrix(0,maxgrid,2)
    gg3[,1]=gg1;gg3[,2]=gg2
    #gg2= seq((gg[1,2]),(gg[2,2]),length=maxgrid)
    gg=gg3
    #print(dim(gg))
    #	Lambda_min<- rat*big_lambda

    #lambda_i<- exp(seq(log(big_lambda),log(big_lambda*rat),length=maxgrid))

    #lambda_i[1]<-big_lambda;lambda_i[maxgrid]<-Lambda_min


    #print(lambda_i)

  }else{
    lambda_i<-my_lambda
    gg1= gg#exp(seq((gg[1]),(gg[2]),length=nlambda))
    #gg2= exp(seq((gg[1,2]),(gg[2,2]),length=nlambda))
    gg=gg1
  }

  lam_list = list()
  beta_0_list = list()
  theta_0_list = list()
  beta_list = list()
  theta_list = list()
  obj=c()
  n_main_terms=c()
  non_zero_theta=c()
  my_obj=list()

  # y.orig=y
  #y=y.orig/sqrt(dot(y.orig,y.orig))
  # my_W_hat<-matrix(0,N,(p*(K+1)))
  # for (i in 1:N) {
  #   my_W_hat[i,c(1:p)]<-W_hat[[i]][,1]
  #   my_W_hat[i,c((p+1):(p*(K+1)))]<-t(W_hat[[i]][,c(2:(K+1))])
  # }

  my_W_hat<-generate.my.w(X=X,Z=Z, quad = TRUE)

  svd.w<- svd(my_W_hat)
  svd.w$tu<- t(svd.w$u)
  svd.w$tv<- t(svd.w$v)


  rho1=rho









  #print(c(length(lam_list),length(n_main_terms),length(non_zero_theta),length(obj)  ) )

  #  for (h in 1:nlambda) {
  #    res_dual <- 0    # dual residual
  #    res_pri <- 0    # primal residual
  #
  #  lambda=lambda_i[h]
  #
  #  #print(lambda)
  #  lambda=lambda*N
  #  beta0 = 0.0#estimates$Beta0
  #  theta0 = matrix(0,K)
  #  beta = matrix(0,p)
  #  V=matrix(0,p,2*(1+K))
  #  O=matrix(0,p,2*(1+K))
  #  Big_beta<-matrix(0,p,2*(1+K))
  #  beta_hat<-matrix(0,p,(1+K))
  #  G=matrix(0,2*(1+K),(1+K))
  #  diag_G<-matrix(0,(K+1),(K+1))
  #  diag(diag_G)<-1
  #  for (i in 1:(K+1)) {
  #    G[i,]<-diag_G[i,]
  #  }
  #  for (i in (K+3):(2*(K+1))) {
  #    G[i,]<-diag_G[(i-(K+1)),]
  #  }
  #
  #
  #  theta = matrix(0,p,K)
  #  Q=matrix(0,p,(1+K))
  #  P=matrix(0,p,(1+K))
  #  my_G=matrix(0,2*(K+1),(p+p*K))
  # my_G[,c(1:p)]=G[,1]
  # my_G[,c((p+1):(p+p*K))]=G[,c(2:(K+1))]
  # new_G<-diag(t(my_G)%*%my_G)
  # V_old<-V;Q_old<-Q
  # my_values<-main_function(rho1,max_it,W,W_hat,my_W_hat,y,N,p,K,e.abs, e.rel,alpha,lambda,alph,my_print=my_print)### iteration
  # beta=unlist(my_values[1]);theta=matrix(unlist(my_values[2]),p,K);converge=unlist(my_values[3])
  #   if(converge==F){
  #       break
  #   }
  # next(h)}## first iteration
  #
  #
  ###############################################################################
  # if(converge==F | h==nlambda ){
  #     if(h==nlambda){h<-h+1}
  #     if(h>1){
  # my_G=matrix(0,2*(K+1),(p+p*K))
  # my_G[,c(1:p)]=G[,1]
  # my_G[,-c(1:p)]=G[,c(2:(K+1))]
  # new_G<-diag(t(my_G)%*%my_G)
  TT<-tree

  C<-TT$Tree
  CW<-TT$Tw


  D=dim(y)[2]


  ### for response groups ###############################################################

  #
  input<-1:(dim(y)[2]*nrow(C))
  multiple_of_D = (input %% dim(y)[2]) == 0

  I<-matrix(0,nrow = nrow(C)*dim(y)[2],ncol = dim(y)[2])

  II<-input[multiple_of_D]
  diag(I[c(1:dim(y)[2] ),])<-C[1,]*(CW[1])

  c_count<-2
  for (e in II[-length(II)]) {
    #	print(e+1)
    #print(diag(I[c((e+1):( c_count*dim(y)[2]) ),]))
    diag(I[c((e+1):( c_count*dim(y)[2]) ),]) <-C[c_count,]*(CW[c_count])
    c_count= 1+c_count
  }
  #
  # #DD=diag((new_G+1),p+p*(K)) #### denominator part for GG=1
  #new_II<-matrix(0,(p+p*K), D) ### number of groups for each response
  new_I=diag(t(I)%*%I)
  #print(I)
  # for (ss in 1:D) {
  #   new_II[,ss]<-new_I[ss]
  #
  # }
  # #
  #
  new_G<-matrix(0,(p+p*K))
  new_G[c(1:p)]<-1;new_G[-c(1:p)]<-2
  new_G[c(1:p)]<-rho*(1+new_G[c(1:p)]);new_G[-c(1:p)]<-rho*(1+new_G[-c(1:p)])

  # #DD=Diagonal(x=new_G) #### denominator part for GG=1
  #
  # DD= Diagonal(x=new_G)#### denominator part for GG=1
  # #DD=matrix(DD, nrow=(p+p*K),ncol = (p+p*K))
  #
  #
  # #invmat <- solve(t(my_W_hat)%*%(my_W_hat)/N+ rho1*DD)
  #
  #
  # XtY <- crossprod(my_W_hat, y)
  #
  #
  #

  invmat<-list() #denominator of the beta estimates
  for (rr in 1:D) {

    DD1<-rho1*(new_I[rr]+1)

    #DD1<-rho1*(obv_beta_matrix)

    DD2<-new_G+DD1

    #DD2[c(1:p)]<-DD2[c(1:p)]+DD1






    #print(rbind(part4[3,]-as.matrix(part4[3,]) ))
    invmat[[rr]] <-DD2# Matrix::chol2inv( Matrix::chol(new_sparse) )
    #print(dim(invmat))
    # # #Matrix::chol2inv()


  }









  beta0 = matrix(0,1,D)#estimates$Beta0
  theta0 = matrix(0,K,D)
  beta =  (matrix(0,p,D))
  beta_hat<-(matrix(0,p+p*(K),D))
  V=(array(0,c(p,2*(1+K),D) ))
  O=(array(0,c(p,2*(1+K),D) ))
  E<-(matrix(0,dim(y)[2]*nrow(C),(p+p*K))) #response auxiliary
  EE<-(array(0,c(p,(1+K),D) ))



  # auxiliary variables for the L1 norm####

  theta =( array(0,c(p,K,D)))
  Q=(array(0,c(p,(1+K),D) ))
  P=(array(0,c(p,(1+K),D) ))
  H<-(matrix(0,dim(y)[2]*nrow(C),(p+p*K)))  # response multiplier
  HH<-(array(0,c(p,(1+K),D) ))
  if(is.null(my_lambda)){
    lam<-matrix(0,nlambda,dim(y)[2])
    for (i in 1:dim(y)[2] ) {
      lam[,i]<-lambda_i[[i]]
    }
  }else{
    lam=lambda_i
  }





  r_current = y#-model(beta0,theta0,beta=beta_hat, theta, X=W_hat, Z)
  b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
  beta0<-b$beta0
  theta0<-b$theta0

  new_y<- y-( matrix(1,N)%*%beta0+Z%*%((theta0)))

  XtY <- crossprod((my_W_hat), (new_y))


  cl1<-cl
  #registerDoMC(cl1)
  if(parallel){
    #registerDoSEQ()
    cl = makeCluster(cl1,type="SOCK")
    registerDoParallel(cl)

    #on.exit(stopCluster(cl))
    getDoParRegistered()
    #registerDoParaqqqllel(numCores)

    my_values<- foreach (i=1:nlambda,.packages='MADMMplasso', .combine=rbind) %dopar% {
      admm.MADMMplasso(beta0=beta0,theta0=theta0,beta=beta,beta_hat=beta_hat,theta=theta,rho1,X,Z,max_it,W_hat=my_W_hat,XtY,y,N,p,K,e.abs, e.rel,alpha, lambda=lam[i,],alph,svd.w=svd.w,tree = tree,my_print=my_print,invmat=invmat,V=V,Q=Q,E=E,EE=EE,O=O,P=P,H=H,HH=HH,cv=cv,gg=gg[i,])
    }
    stopCluster(cl)

  }else if(parallel==F & pal==0){
    my_values=   lapply(seq_len(nlambda),
                        function(g)( admm.MADMMplasso(beta0=beta0,theta0=theta0,beta=beta,beta_hat=beta_hat,theta=theta,rho1,X,Z,max_it,W_hat=my_W_hat,XtY,y,N,p,K,e.abs, e.rel,alpha, lambda=lam[g,],alph,svd.w=svd.w,tree = tree,my_print=my_print,invmat=invmat,V=V,Q=Q,E=E,EE=EE,O=O,P=P,H=H,HH=HH,cv=cv,gg=gg[g,])      ))

  }

  repeat_loop=0
  hh=1
  while (hh<=nlambda) {
    # print(hh)
    res_dual <- 0    # dual residual
    res_pri <- 0    # primal residual
    lambda=lam[hh,]
    #  print(lambda)
    #print(lambda)

    # if(pal==1){
    #
    #   beta0 = matrix(0,1,D)#estimates$Beta0
    #   theta0 = matrix(0,K,D)
    #   beta = matrix(0,p,D)
    #   beta_hat<-matrix(0,p+p*(K),D)
    #   V=array(0,c(p,2*(1+K),D) )
    #   O=array(0,c(p,2*(1+K),D) )
    #   E<-matrix(0,dim(y)[2]*nrow(C),(p)) #response auxiliary
    #   EE<-matrix(0,(p),dim(y)[2])
    #
    #
    #
    #   # auxiliary variables for the L1 norm####
    #
    #   theta = array(0,c(p,K,D))
    #   Q=array(0,c(p,(1+K),D) )
    #   P=array(0,c(p,(1+K),D) )
    #   H<-matrix(0,dim(y)[2]*nrow(C),(p))  # response multiplier
    #   HH<-matrix(0,(p),dim(y)[2])
    #
    # }
    #


    start_time <- Sys.time()
    if(pal==1){
      my_values<-	admm.MADMMplasso(beta0=beta0,theta0=theta0,beta=beta,beta_hat=beta_hat,theta=theta,rho1,X,Z,max_it,W_hat=my_W_hat,XtY,y,N,p,K,e.abs, e.rel,alpha, lambda=lambda,alph,svd.w=svd.w,tree = tree,my_print=my_print,invmat=invmat,V=V,Q=Q,E=E,EE=EE,O=O,P=P,H=H,HH=HH,cv=cv,gg=gg[hh,])

      beta=my_values$beta;theta=my_values$theta;converge=my_values$converge;my_obj[[hh]]<-list(my_values$obj);beta0=my_values$beta0;theta0=my_values$theta0### iteration
      V=my_values$V;Q=my_values$Q;O=my_values$O;P=my_values$P;E=my_values$E;H=my_values$H;beta_hat=my_values$beta_hat; y_hat<-my_values$y_hat
    }
    cost_time <- Sys.time() - start_time
    print(cost_time)
    #
    if(parallel==T & pal==0){
      beta=my_values[hh,]$beta;theta=my_values[hh,]$theta;converge=my_values[hh,]$converge;my_obj[[hh]]<-list(my_values[hh,]$obj);beta0=my_values[hh,]$beta0;theta0=my_values[hh,]$theta0### iteration
      V=my_values[hh,]$V;Q=my_values[hh,]$Q;O=my_values[hh,]$O;P=my_values[hh,]$P;beta_hat=my_values[hh,]$beta_hat;y_hat<-my_values[hh,]$y_hat
      # beta=my_values[[hh]]$beta;theta=my_values[[hh]]$theta;converge=my_values[[hh]]$converge;my_obj[[hh]]<-list(my_values[[hh]]$obj);beta0=my_values[[hh]]$beta0;theta0=my_values[[hh]]$theta0### iteration
      # V=my_values[[hh]]$V;Q=my_values[[hh]]$Q;O=my_values[[hh]]$O;P=my_values[[hh]]$P;beta_hat=my_values[[hh]]$beta_hat;y_hat<-my_values[[hh]]$y_hat
    }else if(parallel==F & pal==0){
      beta=my_values[[hh]]$beta;theta=my_values[[hh]]$theta;converge=my_values[[hh]]$converge;my_obj[[hh]]<-list(my_values[[hh]]$obj);beta0=my_values[[hh]]$beta0;theta0=my_values[[hh]]$theta0### iteration
      V=my_values[[hh]]$V;Q=my_values[[hh]]$Q;O=my_values[[hh]]$O;P=my_values[[hh]]$P;beta_hat=my_values[[hh]]$beta_hat;y_hat<-my_values[[hh]]$y_hat

    }


    #	my_values<-main_function(rho1,X,Z,max_it,W_hat,y,N,p,K,e.abs, e.rel,alpha,lambda,alph,svd.w=svd.w,tree = tree,my_print=my_print)### iteration

    #function(rho1,max_it,W,W_hat,my_W_hat,y,N,p,K,e.abs, e.rel,alpha,lambda,alph,my_print=TRUE,svd.w=svd.w)

    beta1=as(beta*(abs(beta)>tol),"sparseMatrix");theta1=as.sparse3Darray(theta*(abs(theta)>tol));beta_hat1=as(beta_hat*(abs(beta_hat)>tol),"sparseMatrix")
    #beta1=beta;theta1=theta;beta_hat1=beta_hat



    n_interaction_terms = count_nonzero_a((theta1))

    n_main_terms = (c(n_main_terms,count_nonzero_a((beta1)) ) )







    #print(c(length(y),length(model(beta0, theta0, beta, theta, X, Z))))

    # y_hat<-model(beta0, theta0, beta=beta_hat1, theta1, X=my_W_hat, Z)

    obj1<- (sum( as.vector((y - y_hat)^2) ) )/(D*N)#objective(beta0,theta0,
    #  beta,theta,X,Z,
    #  y,W_hat,alpha,lambda,p,K,N)
    obj<-c(obj, obj1)


    non_zero_theta<- (c(non_zero_theta,n_interaction_terms ))
    # print(lambda)
    lam_list<-(c(lam_list,lambda))
    #print(c(length(lam_list),length(n_main_terms),length(non_zero_theta),length(obj)  ) )



    #beta_0_list<-(c(beta_0_list,beta0))
    #theta_0_list[[(hh+1)]]<-list(matrix(theta0,1,K))
    #beta_list[[(hh+1)]]<-beta
    #theta_list[[(hh+1)]]<-list(as.matrix(theta,p,K))

    BETA0[[hh]]<-beta0
    THETA0[[hh]]<-theta0
    BETA[[hh]]<-as(beta1,"sparseMatrix")
    BETA_hat[[hh]]<-as(beta_hat1,"sparseMatrix")

    Y_HAT[[hh]]<-y_hat
    # print(beta1)
    # print(BETA[[b]][,i])
    THETA[[hh]]<-as.sparse3Darray(theta1)


    if(cv==F){
      if(hh==1){
        print(c(hh,(n_main_terms[hh]),non_zero_theta[hh] , obj1))
      }else{
        print(c(hh,(n_main_terms[hh]),non_zero_theta[hh] ,obj[hh-1], obj1))
      }
    }else{
      if(hh==1){
        print(c(hh,(n_main_terms[hh]),non_zero_theta[hh] , obj1))
      }else{
        print(c(hh,(n_main_terms[hh]),non_zero_theta[hh] ,obj[hh-1], obj1))
      }

    }

    #print(length(obj))
    # print(obj)





    hh=hh+1}### lambda
  #  }
  #}

  remove(invmat)
  remove(V); remove(E);remove(H);remove(Q);remove(P);remove(O)
  remove(my_values);remove(my_W_hat)

  obj[1]<-obj[2]
  # print(c(length(lam_list),length(n_main_terms),length(non_zero_theta),length(obj)  ) )

  pred<-data.frame(Lambda=lam,nzero=n_main_terms,nzero_inter=non_zero_theta,OBJ_main=obj)
  out=list(beta0=BETA0,beta=BETA,BETA_hat=BETA_hat,theta0=THETA0,theta=THETA,path=pred,Lambdas=lam,non_zero=n_main_terms,LOSS=obj,it.obj=my_obj,Y_HAT=Y_HAT,gg=gg)
  class(out)="MADMMplasso"
  # Return results
  return (out)
}


convNd2T <- function(Nd, w, w_max){
  # Nd : node list
  # w : a vector of weights for internal nodes
  # Tree : VxK matrix\
  #	V is the number of leaf nodes and internal nodes
  #	K is the number of tasks
  #	Element (v,k) is set to 1 if task k has a membership to
  #	the cluster represented by node v. Otherwise, it's 0.
  # Tw : V vector

  #===========================
  find_leaves <- function(Nd, ch, K, Jt, w, Tw){

    for(ii in 1:length(ch)){
      if(Nd[ch[ii], 2] > K){
        leaves0 <- find_leaves(Nd, which(Nd[,1] == Nd[ch[ii], 2]), K, Jt, w, Tw)
        Jt <- leaves0$Jt
        Tw <- leaves0$Tw
      }else
        Jt <- c(Jt, Nd[ch[ii], 2])
    }

    Tw[Nd[ch, 2]] <- Tw[Nd[ch, 2]] * w

    return(list(Jt=Jt, Tw=Tw))
  }
  #===========================

  # of leaf nodes
  K <- Nd[1,1] - 1
  #V = Nd(size(Nd,1),1);
  #V = Nd(size(Nd,1),1)-1;		# without the root
  if(sum(w < w_max)<1){
    V <- 1 + K
  }else{
    ind0 <- which(w < w_max)    # only the internal nodes with w<w_max
    V <- ind0[length(ind0)] + K
  }

  # for leaf nodes
  I <- 1:K
  J <- 1:K

  Tw <- rep(1, V)

  # for internal nodes
  for(i in (K+1):V){
    Jt <- NULL

    Tw[i] <- Tw[i] * (1 - w[i-K])
    leaves0 <- find_leaves(Nd, which(Nd[,1] == i), K, Jt, w[i-K], Tw)
    Jt <- leaves0$Jt
    Tw <- leaves0$Tw

    I <- c(I, rep(1,length(Jt)) * i)
    J <- c(J, Jt)
  }

  Tree <- sparseMatrix(i=I, j=J, x=rep(1, length(I)), dims=c(V, K))

  return(list(Tree=Tree, Tw=Tw))

}


convH2T <- function(H, w_max){
  K <- dim(H)[1] + 1
  Nd <- cbind(rep((K+1):(2*K-1), each = 2), as.vector(t(H[,1:2])))
  W_norm <- H[,3]/max(H[,3])
  conv0 <- convNd2T(Nd, W_norm, w_max)

  return(conv0)
}

fastCorr <- function(A){
  C <- crossprod(scale(A))/(dim(A)[1]-1)
  return(C)
}

#' Fit the hierarchical tree structure
#' @param y  N by D matrix of response variables
#' @param h is the tree cut off
#' @return  A trained the tree with the following components:
#' Tree: the tree matrix stored in 1's and 0's
#' Tw: tree weights assocuated with the tree matrix. Each weight corresponds to a row in the tree matrix.
#' @export
tree.parms <- function(y=y, h=.7){
  m <- dim(y)[2]
  myDist0 <- 1 - abs(fastCorr(y))
  myDist <- myDist0[lower.tri(myDist0)]
  a0 <- dist(t(y))
  a0[1:length(a0)] <- myDist
  # hierarchical clustering for multivariate responses
  myCluster_0 <- hclust(a0, method = "complete")
  myCluster <- cbind(ifelse(myCluster_0$merge < 0, - myCluster_0$merge, myCluster_0$merge + m), myCluster_0$height)

  conv0 <- convH2T(myCluster, h)
  Tree <- conv0$Tree
  if(is.null(dim(Tree)))
    Tree <- matrix(Tree, nrow=1)
  Tw <- conv0$Tw
  idx <- c(apply(Tree,1,sum) == 1)
  Tree <- Tree[!idx,]
  if(is.null(dim(Tree)))
    Tree <- matrix(Tree, nrow=1)
  Tw <- Tw[!idx]

  #no_group<-	which(colSums(Tree)==0)
  #if(length(no_group)!=0){
  #  tree_matrix<-Matrix(0,(nrow(Tree)+length(no_group)),dim(y)[2],sparse = T )
   # tree_matrix[c(1:nrow(Tree)),]<-Tree
   # count=nrow(Tree)+1
   # for (i in no_group) {

    #  tree_matrix[count,i]<-1
    #  count=count+1
   # }

  #}else{tree_matrix=Tree}
  #tree_weight<-rep(0,length(Tw)+length(no_group))
  #tree_weight[1:length(Tw)]<-Tw;tree_weight[-c(1:length(Tw))]<-1
  out=list(Tree=Tree, Tw=Tw,h_clust=myCluster_0,y.colnames=colnames(y))

  return(out )
}


#'  Simulate data for the model
#' @param p column for X which is the main effect
#' @param n number of observations
#' @param m number of responses
#' @param nz TODO: fill in (number of modifiers?)
#' @param rho  TODO: fill in (correlation between the main effect and modifiers?)
#' @param B.elem  TODO: fill in (the proportion of non-zero elements in beta?)
#' @return  The simulated data with the following components:
#'  Beta: matrix of actual beta coefficients  p by D
#'  Theta: a D by p by K array of actual theta coefficients
#'  Y: a N by D matrix of response variables
#'  X: a N by p matrix of covariates
#'  Z: a N by K matrix of modifiers
#' @export
sim2 <- function(p=500,n=100,m=24,nz=4,rho=.4,B.elem=0.5){
  b<-10
  if(!is.na(p)){
    # generate covariance matrix
    Ap1<-matrix(rep(rho,(p/b)^2),nrow=p/b)
    diag(Ap1)<-rep(1,p/b)
    # Ap2<-matrix(rep(rho,(p[2]/b)^2),nrow=p[2]/b)
    #diag(Ap2)<-rep(1,p[2]/b)
    #Bp12<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[1]/b)
    #Bp21<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[2]/b)
    Xsigma1<-Ap1
    #Xsigma2<-Ap2
    # Xsigma12<-Bp12
    #Xsigma21<-Bp21
    for(i in 2:b){
      Xsigma1<-bdiag(Xsigma1,Ap1)
      # Xsigma12<-bdiag(Xsigma12,Bp12)
      #Xsigma2<-bdiag(Xsigma2,Ap2)
      #Xsigma21<-bdiag(Xsigma21,Bp21)
    }

    Xsigma<-rbind(cbind(Xsigma1))
    X<-	MASS::mvrnorm(n,mu=rep(0,p),Sigma=Xsigma)
    #X1<-X[,1:p[1]]
    #X2<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    #X[,(p[1]+1):(p[1]+p[2])]<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    # generate uncorrelated error term
    esd<-diag(m)
    e<-MASS::mvrnorm(n,mu=rep(0,m),Sigma=esd)

    ## generate beta1 matrix
    Beta1<-matrix(0,nrow=m,ncol=p)
    theta<-array(0,c(p,nz,m))
    Beta1[,1]<-B.elem

    for(i in 1:2){
      Beta1[((i-1)*m/2+1):(i*m/2),(1+(i-1)*2+1):(1+i*2)]<-B.elem
      #theta[(1+(i-1)*2+1),1,((i-1)*m/2+1):(i*m/2)]<- 0.6
    }
    for(i in 1:4){
      Beta1[((i-1)*m/4+1):(i*m/4),(1+4*2+(i-1)*4+1):(1+4*2+i*4)]<-B.elem
      #theta[(1+4*2+(i-1)*4+1),c(2),((i-1)*m/4+1):(i*m/4)]<- 0.6
    }
    for(i in 1:8){
      Beta1[((i-1)*m/8+1):(i*m/8),(1+2*2+4*4+(i-1)*8+1):(1+2*2+4*4+i*8)]<-B.elem
      #theta[(1+2*2+4*4+(i-1)*8+1),3,((i-1)*m/8+1):(i*m/8)]<- 0.6
    }

    theta[30,1,1]<-0.6;theta[31,2,1]<-0.6;theta[32,3,1]<- -0.6;theta[33,4,1]<- -0.6;
    theta[30,1,2]<-0.6;theta[31,2,2]<-0.6;theta[32,3,2]<- -0.6;theta[33,4,2]<- -0.6;

    theta[35,1,5]<-0.6;theta[36,2,5]<-0.6;theta[37,3,5]<- -0.6;theta[38,4,5]<- -0.6;
    theta[35,1,6]<-0.6;theta[36,2,6]<-0.6;theta[37,3,6]<- -0.6;theta[38,4,6]<- -0.6;

    theta[40,1,8]<-0.6;theta[41,2,8]<-0.6;theta[42,3,8]<- -0.6;theta[43,4,8]<- -0.6;
    theta[40,1,9]<-0.6;theta[41,2,9]<-0.6;theta[42,3,9]<- -0.6;theta[43,4,9]<- -0.6;

    theta[48,1,10]<-0.6;theta[49,2,10]<-0.6;theta[50,3,10]<- -0.6;theta[51,4,10]<- -0.6;
    theta[48,1,11]<-0.6;theta[49,2,11]<-0.6;theta[50,3,11]<- -0.6;theta[51,4,11]<- -0.6;

    theta[57,1,13]<-0.6;theta[58,2,15]<-0.6;theta[59,3,15]<- -0.6;theta[60,4,15]<- -0.6;
    theta[57,1,15]<-0.6;theta[58,2,17]<-0.6;theta[59,3,17]<- -0.6;theta[60,4,17]<- -0.6;

    theta[63,1,16]<-0.6;theta[64,2,16]<-0.6;theta[65,3,16]<- -0.6;theta[66,4,16]<- -0.6;
    theta[63,1,17]<-0.6;theta[64,2,17]<-0.6;theta[65,3,17]<- -0.6;theta[66,4,17]<- -0.6;

    theta[80,1,21]<-0.6;theta[81,2,21]<-0.6;theta[82,3,21]<- -0.6;theta[83,4,21]<- -0.6;
    theta[80,1,22]<-0.6;theta[81,2,22]<-0.6;theta[82,3,22]<- -0.6;theta[83,4,22]<- -0.6

    ## generate beta2 matrix
    #   Beta2<-matrix(0,nrow=m,ncol=p[2])
    #   Beta2[,1]<-B.elem[2]
    #   for(i in 1:3){
    #     Beta2[((i-1)*m/3+1):(i*m/3),(1+(i-1)*2+1):(1+i*2)]<-B.elem[2]
    #   }
    #   for(i in 1:6){
    #     Beta2[((i-1)*m/6+1):(i*m/6),(1+2*3+(i-1)*4+1):(1+2*3+i*4)]<-B.elem[2]
    #   }
    #   for(i in 1:12){
    #     Beta2[((i-1)*m/12+1):(i*m/12),(1+2*3+4*6+(i-1)*8+1):(1+2*3+4*6+i*8)]<-B.elem[2]
    #   }
    #   Beta<-t(cbind(Beta1, Beta2))
    # }else{
    #   cat("Ooops!!! Please specify 2-dim p vector, for example p=c(500,150)\n")
  }
  Beta<-t((Beta1))

  mx=colMeans(X)

  sx=sqrt(apply(X,2,var))
  X=scale(X,mx,sx)
  Z= matrix(rbinom(n = n*nz, size = 1, prob = 0.5), nrow = n, ncol = nz)
  # Z =matrix(rnorm(n*nz),n,nz)
  #mz=colMeans(Z)
  # sz=sqrt(apply(Z,2,var))

  #Z=scale(Z,mz,sz)

  pliable = matrix(0,n,m)
  for (ee in 1:m) {
    pliable[,ee]<-	compute_pliable(X, Z, theta[,,ee])

  }
  # meta_matrix<-matrix(0,p[1]+p[1]*nz,m)
  # for (i in 1:m) {
  #   meta_matrix[c(1:p),i]<-Beta[,i]
  #   meta_matrix[-c(1:p),i]<-as.vector(theta[,,i] )
  # }


  #X=matrix(as.numeric(X),N,p)
  # X= scale(X)
  Y<-X%*%Beta+pliable+e
  #Y<-X%*%Beta+e

  # Z <- matrix(rbinom(n = n*nz, size = 1, prob = 0.5), nrow = n, ncol = nz)
  # Y=Y+rep(rowSums(cbind(0.6*X[,1]*Z[,1],0.6*X[,3]*Z[,2],0.6*X[,10]*Z[,3],0.6*X[,12]*Z[,4])),m)
  out=list(Y=Y, X=X,Z=Z, Beta=Beta,Theta=theta, e=e, p=p)
  #class(out)="sim2"
  return(out)
}






# Rcpp::cppFunction(
#   depends = "RcppArmadillo",
#   code    = "
#     arma::mat tcrossprod_cpp(const arma::mat& x, const arma::mat& y) {
#       return(x * y.t());
#     }
#   "
# )



errfun.gaussian<-function(y,yhat,w=rep(1,length(y))){  ( w*(y-yhat)^2) }


#' Carries out cross-validation for  a  pliable lasso model over a path of regularization values
#' @param fit  object returned by the pliable function
#' @param X  N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may
#' represent quantitative or categorical variables, or a mixture of the two.
#'  Categorical varables should be coded by 0-1 dummy variables: for a k-level
#' variable, one can use either k or k-1  dummy variables.
#' @param y N by D-matrix of responses. The X and Z variables are centered in
#' the function. We recommmend that x and z also be standardized before the call
#' @param nfolds  number of cross-validation folds
#' @param foldid  vector with values in 1:K, indicating folds for K-fold CV. Default NULL
#' @param alpha TODO: add parameter description
#' @param lambda TODO: add parameter description
#' @param max_it TODO: add parameter description
#' @param e.abs TODO: add parameter description
#' @param e.rel TODO: add parameter description
#' @param nlambda TODO: add parameter description
#' @param rho TODO: add parameter description
#' @param my_print TODO: add parameter description
#' @param alph TODO: add parameter description
#' @param foldid TODO: add parameter description
#' @param parallel TODO: add parameter description
#' @param pal TODO: add parameter description
#' @param gg TODO: add parameter description
#' @param TT TODO: add parameter description
#' @param tol TODO: add parameter description
#' @param cl TODO: add parameter description
#' @return  results containing the CV values
#' @examples
#' \dontrun{
#' # Train the model
#' # generate some data
#' set.seed(1235)
#' N = 100 ; p =50;nz=4; K=nz
#' X <- matrix(rnorm(n = N * p), nrow = N, ncol = p)
#' mx=colMeans(X)
#' sx=sqrt(apply(X,2,var))
#' X=scale(X,mx,sx)
#' X=matrix(as.numeric(X),N,p)
#' Z =matrix(rnorm(N*nz),N,nz)
#' mz=colMeans(Z)
#' sz=sqrt(apply(Z,2,var))
#' Z=scale(Z,mz,sz)
#' beta_1 <- rep(x = 0, times = p)
#' beta_2<-rep(x = 0, times = p)
#' beta_3<-rep(x = 0, times = p)
#' beta_4<-rep(x = 0, times = p)
#' beta_5<-rep(x = 0, times = p)
#' beta_6<-rep(x = 0, times = p)
#'
#'
#' beta_1[1:5] <- c(2, 2, 2, 2,2)
#' beta_2[1:5]<-c(2, 2, 2, 2,2)
#' beta_3[6:10]<-c(2, 2, 2, -2,-2)
#' beta_4[6:10] <- c(2, 2, 2, -2,-2)
#' beta_5[11:15] <- c(-2,  -2,-2, -2,-2)
#' beta_6[11:15] <- c(-2, -2, -2, -2,-2)
#'
#' Beta<-cbind(beta_1,beta_2,beta_3,beta_4,beta_5,beta_6)
#' colnames(Beta)<-c(1:6)
#'
#'
#' theta<-array(0,c(p,K,6))
#' theta[1,1,1]<-2;theta[3,2,1]<-2;theta[4,3,1]<- -2;theta[5,4,1]<- -2;
#' theta[1,1,2]<-2;theta[3,2,2]<-2;theta[4,3,2]<- -2;theta[5,4,2]<- -2;
#' theta[6,1,3]<-2;theta[8,2,3]<-2;theta[9,3,3]<- -2;theta[10,4,3]<- -2;
#' theta[6,1,4]<-2;theta[8,2,4]<-2;theta[9,3,4]<- -2;theta[10,4,4]<- -2;
#' theta[11,1,5]<-2;theta[13,2,5]<-2;theta[14,3,5]<- -2;theta[15,4,5]<- -2;
#' theta[11,1,6]<-2;theta[13,2,6]<-2;theta[14,3,6]<- -2;theta[15,4,6]<- -2
#'

#' pliable = matrix(0,N,6)
#' for (e in 1:6) {
#' pliable[,e]<-	compute_pliable(X, Z, theta[,,e])
#' }
#'
#' esd<-diag(6)
#' e<-MASS::mvrnorm(N,mu=rep(0,6),Sigma=esd)
#' y_train<-X%*%Beta+pliable+e
#' y=y_train
#' #colnames(y)<-c(1:6)
#' colnames(y)<- c( paste("y",1:(ncol(y)),sep = "") )
#' TT=tree.parms(y)
#' plot(TT$h_clust)
#' gg1=matrix(0,2,2)
#' gg1[1,]<-c(0.02,0.02)
#' gg1[2,]<-c(0.02,0.02)
#' nlambda = 50
#' e.abs=1E-4
#' e.rel=1E-2
#' alpha=.2
#' tol=1E-3
#' fit <- MADMMplasso(
#'   X, Z, y, alpha=alpha, my_lambda=NULL, lambda_min=0.001, max_it=5000,
#'   e.abs=e.abs, e.rel=e.rel, maxgrid=50, nlambda=nlambda, rho=5,tree=TT,
#'   my_print=FALSE, alph=1, parallel=FALSE, pal=1, gg=gg1, tol=tol, cl=6
#' )
#' gg1=fit$gg
#'
#' cv_admp <- cv.MADMMplasso(
#'   fit, nfolds=5, X, Z, y, alpha=alpha, lambda=fit$Lambdas, max_it=5000,
#'   e.abs=e.abs, e.rel=e.rel, nlambda, rho=5, my_print=FALSE, alph=1,
#'   foldid=NULL, parallel=FALSE, pal=1, gg=gg1, TT=TT, tol=tol
#' )
#' plot(cv_admp)
#' }
#' @export
cv.MADMMplasso<-function(fit,nfolds,X,Z,y,alpha=0.5,lambda=fit$Lambdas,max_it=50000,e.abs=1E-3,e.rel=1E-3,nlambda, rho=5,my_print=FALSE,alph=1,foldid=NULL,parallel=TRUE,pal=0,gg=c(7,0.5),TT,tol=1E-4,cl=2){
  BIG=10e9
  no<-nrow(X)
  ni<-ncol(X)
  nz<-ncol(Z)
  ggg=vector("list",nfolds)

  yhat=array(NA,c(no,dim(y)[2],length(lambda[,1]) ))

  # yhat=matrix(0,nfolds,length(result$Lambdas))
  my_nzero<-matrix(0,nfolds,length(lambda[,1]))


  if(is.null(foldid)) foldid = sample(rep(1:nfolds, ceiling(no/nfolds)), no, replace=FALSE)  #foldid = sample(rep(seq(nfolds), length = no))

  nfolds=length(table(foldid))

  status.in=NULL

  for(ii in 1:nfolds){
    print(c("fold,", ii))
    oo=foldid==ii



    ggg[[ii]]<-   MADMMplasso(X=X[!oo,,drop=F],Z=Z[!oo,,drop=F],y=y[!oo,,drop=F],alpha = alpha,my_lambda=lambda,lambda_min=.01,max_it=max_it,e.abs=e.abs,e.rel=e.rel,nlambda=length(lambda[,1]), rho=rho,tree = TT,my_print = my_print,alph=alph,cv=TRUE,parallel = parallel,pal=pal,gg=gg,tol=tol,cl=cl)



    cv_p<-predict.MADMMplasso(ggg[[ii]] ,X=X[oo,,drop=F],Z=Z[oo,],y=y[oo,])
    #PADMM_predict_lasso<-function(object ,X,Z,y,lambda=NULL)

    #  Coeff<-lapply(seq_len(max(y)),
    #               function(j)(matrix(0,ncol(X),nlambda)))

    ggg[[ii]]<-0
    #for (i in 1:nlambda) {
    #  f<-Matrix(unlist(ggg[[ii]]$beta[i]),ncol(X),max(y),sparse = T)
    # for (j in 1:max(y)) {
    #  Coeff[[j]][,i]<-f[,j]
    #}

    # }

    #nnn_zero<-matrix(0,max(y),nlambda)

    #for (i in 1:max(y)) {
    #  q<-matrix(unlist(Coeff[[i]]),ncol(X),nlambda)
    # nnn_zero[i,]<-colSums(q!=0)


    # }
    #non_zero<-matrix(0,nlambda)
    #for (i in 1:nlambda) {
    # non_zero[i]<-max(nnn_zero[,i])
    #}

    #print(cv_p$deviance)
    # my_nzero[ii,]<-non_zero
    # print(unlist(cv_p$y_hat))
    # print(yhat[oo,])
    yhat[oo, , seq(nlambda)] = cv_p$y_hat[, , seq(nlambda)]


    #(result ,X,Z,y,lambda=NULL)
  }

  #print(yhat)
  ym=array(y,dim(yhat))
  #print(ym)
  #err<-matrix(0,length(lambda[,1]),dim(y)[2])
  # for (ii in 1:dim(y)[2]) {
  err<-apply((ym - yhat)^2,c( 1,3), sum)

  #
  mat=err
  outmat = matrix(NA, nfolds, ncol(mat))
  good = matrix(0, nfolds, ncol(mat))
  mat[is.infinite(mat)] = NA
  for (i in seq(nfolds)) {
    mati = mat[foldid == i, , drop = FALSE]
    # wi = weights[foldid == i]
    outmat[i, ] = apply(mati, 2, mean, na.rm = TRUE)
    # good[i, seq(nlams[i])] = 1
  }
  #
  #
  #err= yhat

  #err=outmat

  #err <- lapply(seq_len(length(lambda[,1])),
  #                   function(j) (norm(ym[,,j]-yhat[,,j],type = "F"))^2/(2*N) )

  #err<-matrix(unlist(err),1)
  #non_zero<-matrix(0,nlambda)

  non_zero<-c(fit$path$nzero)

  cvm=(apply(err,2,mean,na.rm=T))/dim(y)[2]
  nn=apply(!is.na(err),2,sum,na.rm=T)
  cvsd=sqrt(apply(err,2,var,na.rm=T)/(dim(y)[2]*nn))

  #cvm<-colMeans(cvm); cvsd<-colMeans(cvsd)
  cvm.nz=cvm
  cvm.nz[non_zero==0]=BIG
  imin=which.min(cvm.nz)
  imin.1se=which(cvm< cvm[imin]+cvsd[imin])[1]

  out=list(lambda=fit$Lambdas,cvm=cvm,cvsd=cvsd,cvup = cvm +
             cvsd, cvlo = cvm - cvsd, nz=c(fit$path$nzero),lambda.min=fit$Lambdas[imin,1],lambda.1se=fit$Lambdas[imin.1se,1])
  class(out)="cv.MADMMplasso"

  return(out)
}



error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
 # range(upper, lower)
}





S_func <- function(x, a) {  # Soft Thresholding Operator
  return(pmax(abs(x) - a,0) * sign(x))
}

#' @title TODO: fill this field
#' @description TODO: fill this field
#' @param X TODO: fill this field
#' @param Z TODO: fill this field
#' @param theta TODO: fill this field
#' @return TODO: fill this field
#' @export
compute_pliable<-function(X, Z, theta){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)

  xz_theta <- lapply(seq_len(p),
                     function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j])
  xz_term<- (Reduce(f = '+', x = xz_theta))


  return(xz_term)


}

model_p<-function(beta0, theta0, beta, theta, X, Z){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  D=dim(beta0)[2]
  #The pliable lasso model described in the paper
  #y ~ f(X)

  #formulated as

  #y ~ b_0 + Z theta_0 + X b + \sum( w_j theta_ji )




  #beta0<-array(0,c(1,1,D))

  #print(D)

  intercepts = matrix(1,N)%*%beta0+Z%*%(theta0)
  #intercepts1<-matrix(0,N,D)
  #intercepts1[,]<-intercepts
  #print(intercepts)
  shared_model = X%*%(beta)
  #shared_model1<-as.vector(t(X))%*%as.vector(t(beta))
 # pliable = matrix(0,N,D)
  #for (e in 1:D) {
  #  pliable[,e]<-	compute_pliable(X, Z, theta[,,e])

  #}



  #apply(theta,compute_pliable,X=X,Z=Z,theta=theta)

  return(intercepts+  shared_model )
}



model_intercept<-function( beta0, theta0, beta, theta, X, Z){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  D=dim(beta0)[2]
  #The pliable lasso model described in the paper
  #y ~ f(X)

  #formulated as

  #y ~ b_0 + Z theta_0 + X b + \sum( w_j theta_ji )




  #beta0<-array(0,c(1,1,D))

  #print(D)

  #intercepts = matrix(1,N)%*%beta0+Z%*%(theta0)
  #intercepts1<-matrix(0,N,D)
  #intercepts1[,]<-intercepts
  #print(intercepts)
  shared_model = X%*%(beta)
  #shared_model1<-as.vector(t(X))%*%as.vector(t(beta))
  # pliable = matrix(0,N,D)
  #for (e in 1:D) {
  #  pliable[,e]<-	compute_pliable(X, Z, theta[,,e])

  #}



  #apply(theta,compute_pliable,X=X,Z=Z,theta=theta)

  return(  shared_model )
}

#




objective<-function(beta0,theta0,beta,theta,X,Z,y,alpha,lambda,p,N,IB,W,beta1){
  #print(length(y))

  loss<- ( norm(y-model_p(beta0, theta0, beta=beta1, theta, X=W, Z), type="F")^2 )
  # print(length(matrix(model(beta0, theta0, beta, theta, X, Z),N,1 )))
  #rbind(matrix(y,nrow  =length(y),ncol =1),matrix(model(beta0, theta0, beta, theta, X, Z),nrow=nrow(X),ncol =1 ))
  mse=(1 / (2*N )) * loss
  #beta<-matrix(beta,p,1);theta<-matrix(theta,p,K)
  #mse = nll_p(beta0,theta0,beta,theta,X,Z,y)

  #lambda=lambda


  l_1<-sum(abs(beta))
  pliable_norm<-matrix(0,dim(y)[2])
  for (ee in 1:dim(y)[2]) {

    beta11<-beta[,ee]; theta11<-theta[,,ee]
    norm_1_l=   lapply(seq_len(p),
                       function(g)(  lambda[ee]*(1-alpha)*(norm(matrix(c(beta11[g],theta11[g,])),type = "F") +norm(matrix(c(theta11[g,])),type = "F") ) + lambda[ee]*alpha* sum(abs(theta11[g,]))      ))
    #


    pliable_norm[ee]<- sum(unlist(norm_1_l))

  }
  objective_l <- mse +(1-alpha)*min(lambda/4)*IB+(1-alpha)*min(lambda/4)*l_1+ sum(pliable_norm)
  #pr<-(sign(y_hat_l))


  # prob_d <-matrix((pr))
  #objective_l=apply(apply(prob_d, 2, FUN="!=", y), 2, sum)/length(y)

  return( objective_l)
}



count_nonzero_a<-function(x){
  if (length(dim(x))==3 ) {
    count1<-matrix(0,dim(x)[3])
    for (ww in 1:dim(x)[3] ) {
      n = sum(x[,,ww]!=0)
      # for (i in 1: length(x[,,ww])){
      #   tet<-x[,,ww]
      #   if ( isTRUE(abs(tet[i])>0)==T ){
      #     n=n + 1
      #   }
      # }

      count1[ww]<-n
    }
    n=max(count1)
  } else{
    count1<-matrix(0,dim(x)[2])
    for (ww in 1:dim(x)[2] ) {
      n = sum(x[,ww]!=0)
      # for (i in 1: length(x[,ww])){
      #   tet<-x[,ww]
      #   if (isTRUE(abs(tet[i])>0)==T ){
      #     n=n + 1
      #   }
      # }

      count1[ww]<-n
    }
    n=max(count1)



  }



  return (n)

}

reg<-function(r,Z){
  K=ncol(Z)
  N=nrow(Z)
  #r=rowMeans(r)
  #r=matrix(r,ncol = 1)
  beta01<-matrix(0,1,ncol(r))
  theta01<-matrix(0,ncol(Z),ncol(r))
  for (e in 1:ncol(r)) {


    #my_one<-matrix(1,nrow(Z))
    #my_w=data.frame(Z)
    #my_w<-as.matrix(my_w)
    #my_inv<-pinv((t(my_w)%*%my_w)/N)
    #my_res<-my_inv%*%( (t(my_w)%*%r[,e])/N )
    #new<- lm(r[,e]~1,na.action=na.exclude,singular.ok = TRUE)


    #beta01[e]<-matrix(my_res[(K+1)])
    new1<- lm(r[,e]~Z,singular.ok = TRUE)
    beta01[e]<-matrix(new1$coefficients[1])
    theta01[,e]<- as.vector(new1$coefficients[-1] )
    #theta01[,e]<- matrix(my_res[c(1:K )])

  }
  # print(beta0)
  # print(theta0)
  return(list(beta0=beta01,theta0=theta01))
}





quick.func<- function(xz = c(),xn){
  as.vector(xz[1:xn]%o%xz[-(1:xn)])
}

#' Generate the matrix W as used in Appendix I for use in the function.
#' @param X TODO: fill in description
#' @param Z TODO: fill in description
#' @param quad TODO: fill in description
#' @export
generate.my.w<- function(X=matrix(),Z=matrix(), quad = TRUE){

  p1<- ncol(X)
  p2<- ncol(Z)

  if(quad == FALSE){
    p<- p1
    if(p1!=p2) stop("To remove quadtratic terms p1 must be equal to p2")
    ind<- (1:p)*p + (2*(1:p)+1)
  }

  #Just in case we have only one oberservation? Not sure why I did this
  if(is.vector(X)) p1<- length(X)
  if(is.vector(Z)) p2<- length(Z)

  #Add the intercept
  x<- X
  z<- cbind(1,Z)

  W<- t(apply(cbind(x,z),1,quick.func,xn= p1))
  if(quad == FALSE){
    W<- W[,-ind]
  }

  return(W)

}


twonorm <- function(x) sqrt(sum(x^2,na.rm = TRUE))




#' Compute predicted values from a fitted pliable  object

#'  Make predictions from a fitted pliable lasso model
#' @param object object returned from a call to pliable
#' @param X  N by p matrix of predictors
#' @param Z  n by nz matrix of modifying variables. These may be observed or
#' the predictions from a supervised learning algorithm that predicts z from
#' test features x  and possibly other features.
#' @param y N by D matrix  of responses.
#' @param lambda  TODO: fill in description
#' @param ... additional arguments to the generic \code{predict()} method
#'  @return  predicted values
#' @export
predict.MADMMplasso<-function(object ,X,Z,y,lambda=NULL, ...){
  lambda.arg=lambda
  if(is.null(lambda.arg))
  { lambda=object$Lambdas[,1]
  isel=1:length(lambda)}

  if(!is.null(lambda.arg)){

    isel=as.numeric(knn1(matrix(object$Lambdas[,1],ncol=1),matrix(lambda.arg,ncol=1),1:length(object$Lambdas[,1])))

  }

  # print(c(isel,length(isel)))

  N <- nrow(X)

  p <- ncol(X)
  #print(Z)
  K <- ncol(as.matrix(Z))
  D=dim(y)[2]
  #print(K)
  #Z=matrix(as.numeric(Z),N,K)
  my_W_hat<-generate.my.w(X=X,Z=Z, quad = TRUE)

  yh=array(0,c(N,D,length(isel)))
  DEV=matrix(NA,length(isel))
  my_theta<-array(0,c(ncol(X),ncol(Z), length(isel)))


  #pBETA0<-matrix(0,nrow = length(isel)); pBETA<-matrix(0,nrow = length(isel)); pTHETA0<-matrix(0,nrow = length(isel)); pTHETA<-matrix(0,nrow =length(isel))


  # pBETA0<-lapply(seq_len(1),
  #                function(j)(matrix(0,nrow = length(isel))))
  #
  # pBETA<-lapply(seq_len(1),
  #               function(j)(matrix(0,nrow=p,ncol=length(isel))))
  #
  # pTHETA0<-lapply(seq_len(1),
  #                 function(j)(matrix(0,nrow=K,ncol=length(isel))))
  #
  # pTHETA<-lapply(seq_len(1),
  #                function(j)(array(0,c(p,K,length(isel)))))
  #
  #
  pBETA0<-lapply(seq_len(length(isel)),
                 function(j)(matrix(0,nrow = (D))))
  #
  pBETA<-lapply(seq_len(length(isel)),
                function(j)(matrix(0,nrow=p,ncol=(D) )) )
  pBETA_hat<-lapply(seq_len(length(isel)),
                    function(j)(matrix(0,nrow=p+p*K,ncol=(D) )) )


  pTHETA0<-lapply(seq_len(length(isel)),
                  function(j)(matrix(0,nrow=K,ncol=(D)  )))

  pTHETA<-lapply(seq_len(length(isel)),
                 function(j)(array(0,c(p,K,(D)  ))))

  pY_HAT<-lapply(seq_len(length(isel)),
                 function(j)( matrix(0,nrow=N,ncol=(D) )  )  )



  ii=0
  for(m in isel){
    ii=ii+1

    # pred_beta<-lapply(seq_len(max(y)),
    #                  function(j)(matrix(0,nrow = 1,ncol = ncol(X))))
    #pred_theta<-lapply(seq_len(max(y)),
    #                    function(j)(matrix(0,nrow = p,ncol = K)))
    #  pred_theta0<-lapply(seq_len(max(y)),
    #                    function(j)(matrix(0,nrow = 1,ncol = K)))

    #  pred_beta0<-lapply(seq_len(max(y)),
    #                  function(j)(0))


    z=m
    #BETA0<-matrix(unlist(object$beta0[z]),1,max(y))


    # BETA<-as.data.frame(object$beta[z])
    # THETA0<-as.data.frame(object$theta0[z])
    # THETA<-object$theta[z]




    n_i<-lapply(seq_len(max(1)),
                function(j)(matrix(0,nrow = N)))
    pr<-lapply(seq_len(max(1)),
               function(j)(matrix(0,nrow = N)))



    # for (x in 1:1) {
    #theta<- matrix(unlist(THETA[[1]][x]),p,K)


    beta0<-object$beta0[[z]]
    beta <- object$beta[[z]]
    beta_hat<-object$BETA_hat[[z]]
    theta <- object$theta[[z]]


    theta0 <- object$theta0[[z]]


    pBETA0[[ii]] <-beta0
    pBETA[[ii]]<-beta
    pBETA_hat[[ii]]<-beta_hat
    pTHETA[[ii]] <-theta
    pTHETA0[[ii]]<-theta0







    # pred_theta[x]<-list(as.matrix(theta,p,K))
    # pred_beta0[x]<-beta0
    # pred_theta0[x]<-theta0
    # pred_beta[x]<-beta
    # beta=matrix(unlist(BETA[x]),p)
    #n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
    n_i<-(model_p(beta0, theta0, beta=beta_hat, theta, X=my_W_hat, Z) )

    #  }



    # pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))

    #v1_d<-1*(v1==d)



    #pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))


    #v1_d<-1*(v1==d)

    # n_i_d<-matrix(unlist(n_i[d]),N)

    # deviance_y1[d]<-y_d*n_i_d[l]
    # deviance_y2[d]<-exp(n_i_d[l])

    # }

    Dev<-( sum( as.vector( (y - n_i)^2 ) ) )/(D*N)


    #Dev1[l]<-  sum(deviance_y2)


    #(y,yhat,w=rep(1,length(y)))
    #}
    #DEV1[i]<-((2)*sum(Dev1))/length(y)
    DEV[ii]<-(Dev)

    # DEV[ii]<- ((2)*sum(Dev))#+
    yh[,,ii]<- as.matrix(n_i)

  }
  out=list(y_hat=yh,beta0=pBETA0,beta=pBETA,beta_hat=pBETA_hat,theta0=pTHETA0,theta=pTHETA,deviance=DEV)
  #class(out)="predict"
  return(out)
}

#' @export
plot.MADMMplasso=
  function(x,...){
    fit=x
    beta<-fit$beta; theta<-fit$theta
    p=nrow(fit$beta[[1]]);K<-nrow(fit$theta0[[1]]);D=ncol(fit$theta0[[1]]);nlambda<-length(fit$Lambdas[,1])

    plotCoeff(beta,theta,error = fit$path$OBJ_main,nz=fit$non_zero,p,K,D,nlambda,Lambda = fit$Lambdas)

  }


plotCoeff=
  function(beta,theta,error,nz,p,K,D,nlambda,Lambda){
    if (nlambda == 1) {
      stop("nlambda must be greater than 1 to produce plot")
    }

    gg=nlambda
    my_beta1<-array(NA,c(p,nlambda,D) )
    for (r in 1:nlambda) {
      # q<-matrix(unlist(Coeff[[i]]),ncol(X),nlambda)




      #result$beta[[r]][1,]
      for (i in 1:D) {


        my_beta1[,r,i]<-beta[[r]][,i]


      }


    }



    for (ii in 1:D) {

      my_beta <- my_beta1[, , ii]

      b <- apply(abs(my_beta), 2, sum)
      b=log(unlist(Lambda[,ii]))
      n=dim(my_beta)[2]
      matplot(b, t(my_beta),type="n",col="red", ylim = range(my_beta),  xlab="Log Lambda", ylab=  ( paste(  "coefficient",ii)))
      axis(side = 3, at =  (as.matrix(b)), labels = paste(as.matrix(nz[c(1:gg)])),
           tick = FALSE, line = 0)


      for (i in 1:(p)) {
        lines(b,(my_beta[i,]),col=i+1, lty=1  )

        text( (min(b-.1)),my_beta[i,n], labels = i,cex=.7)
      }


      my_beta<-(my_beta)
      act=which(rowSums(abs(my_beta))>0)
      theta1<-array(0,c(p,K,(nlambda) ))
      for (i in 1:(nlambda)) {
        theta1[,,i]<-matrix(unlist(theta[[i]][,,ii]),p,K)
      }


      ntheta=  apply( abs(theta1)>1E-3,c(1,3),sum)
      index = b#c(1:gg)#apply(abs(my_beta), 2, sum)+apply(abs(theta1),3,sum)
      sbeta=(my_beta)
      for(j in act){
        for(i in 1:length(index)){

          if(ntheta[j,i]>0) text(index[i],sbeta[j,i],label="x",cex=.7)
        }}


    }
  }



#' @export
plot.cv.MADMMplasso=
  function(x, ...){

    cvobj=x
    xlab = "log(Lambda)"
    plot.args = list(x = log(as.matrix(cvobj$lambda[,1])), y = as.matrix(cvobj$cvm),
                     ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = "Error",
                     type = "n")


    new.args = list(...)
    if (length(new.args))
      plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    error.bars( log(cvobj$lambda[,1]), cvobj$cvup,cvobj$cvlo,
                width = 0.01, col = "darkgrey")
    points(log(as.matrix(cvobj$lambda[,1])), as.matrix(cvobj$cvm), pch = 20,
           col = "red")
    axis(side = 3, at =  log(as.matrix(cvobj$lambda[,1])), labels = paste(as.matrix(cvobj$nz)),
         tick = FALSE, line = 0)

    abline(v= log(cvobj$lambda.min), lty = 3)
    abline(v= log(cvobj$lambda.1se), lty = 3)
    invisible()




  }
