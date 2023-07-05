##### all these functions are called in that particular function#######


S_func <- function(x, a) {  # Soft Thresholding Operator
  return(pmax(abs(x) - a,0) * sign(x))
}

compute_pliable<-function(X, Z, theta){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)

  xz_theta <- lapply(seq_len(p),
                     function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j])
  xz_term<- (Reduce(f = '+', x = xz_theta))

  return(xz_term)


}



model<-function(beta0, theta0, beta, theta, X, Z){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  D=dim(beta0)[2]
  #The pliable lasso model described in the paper
  #y ~ f(X)

  #formulated as

  #y ~ b_0 + Z theta_0 + X b + \sum( w_j theta_ji )




  #beta0<-array(0,c(1,1,D))



  intercepts = (matrix(1,N))%*%beta0+Z%*%((theta0))
  shared_model = X%*%(beta)
  #shared_model1<-as.vector(t(X))%*%as.vector(t(beta))
  pliable = matrix(0,N,D)
  #for (e in 1:D) {
  #  pliable[,e]<-	compute_pliable(X, Z, theta[,,e])

  #}



  #apply(theta,compute_pliable,X=X,Z=Z,theta=theta)

  return( intercepts + shared_model )
}

#

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




library(pracma)

reg<-function(r,Z){
  K=ncol(Z)

  beta01<-matrix(0,1,ncol(r))
  theta01<-matrix(0,ncol(Z),ncol(r))
  for (e in 1:ncol(r)) {


    my_one<-matrix(1,nrow(Z))
    my_w=data.frame(Z,my_one)
    my_w<-as.matrix(my_w)
    my_inv<-pinv(t(my_w)%*%my_w)
    my_res<-my_inv%*%(t(my_w)%*%r[,e])
    # new<- lm(r~1,na.action=na.exclude)
    beta01[e]<-matrix(my_res[(K+1)])
    # new1<- lm(r~Z,singular.ok = TRUE)

    theta01[,e]<- matrix(my_res[c(1:(K))])

  }
  return(list(beta0=beta01,theta0=theta01))
}

twonorm <- function(x) sqrt(sum(x^2,na.rm = TRUE))


############ below are used to generate the data






quick.func<- function(xz = c(),xn){
  as.vector(xz[1:xn]%o%xz[-(1:xn)])
}

#Generate the matrix \widetilde{W} as used in Appendix I for use in the function.
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




convNd2T <- function(Nd, w, w_max){
  # Nd : node list
  # w : a vector of weights for internal nodes
  # Tree : VxK matrix
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





#### I use thye below codes to generate the data##############

library(stats)
library(MASS)
library(Matrix)

set.seed(1235)
N = 100 ; p =50;nz=4; K=nz
X <- matrix(rnorm(n = N * p), nrow = N, ncol = p)

mx=colMeans(X)

sx=sqrt(apply(X,2,var))
X=scale(X,mx,sx)
X=matrix(as.numeric(X),N,p)

#Z <- matrix(rbinom(n = N * K, size = 1, prob = 0.5), nrow = N, ncol = K)
Z =matrix(rnorm(N*nz),N,nz)
mz=colMeans(Z)
sz=sqrt(apply(Z,2,var))

Z=scale(Z,mz,sz)
e=matrix(1,N)
#X <- matrix(rnorm(n = N * p, mean = 0, sd = 1), nrow = N, ncol = p)
#Z <- matrix(rbinom(n = N * K, size = 1, prob = 0.5), nrow = N, ncol = K)

e=matrix(1,N)
beta_1 <- rep(x = 0, times = p); beta_2<-rep(x = 0, times = p);beta_3<-rep(x = 0, times = p);beta_4<-rep(x = 0, times = p);beta_5<-rep(x = 0, times = p);beta_6<-rep(x = 0, times = p)
#beta_1[1:5] <- c(-2, -2, 2, 2,2); beta_2[1:5]<-c(-2, -2, 2, 2,2); beta_3[6:10]<-c(-2, -2, 2, -2,-2); beta_4[6:10] <- c(-2, -2, 2, -2,-2); beta_5[11:15] <- c(-2,  -2,-2, -2,-2);beta_6[11:15] <- c(-2, -2, -2, -2,-2)

beta_1[1:5] <- c(2, 2, 2, 2,2); beta_2[1:5]<-c(2, 2, 2, 2,2); beta_3[6:10]<-c(2, 2, 2, -2,-2); beta_4[6:10] <- c(2, 2, 2, -2,-2); beta_5[11:15] <- c(-2,  -2,-2, -2,-2);beta_6[11:15] <- c(-2, -2, -2, -2,-2)


Beta<-cbind(beta_1,beta_2,beta_3,beta_4,beta_5,beta_6)
colnames(Beta)<-c(1:6)
#beta_1[1:5]<-c(4,-2,2,2,-5)
# coeffs1 <- cbind(beta_1[1]+2*Z[,1], beta_1[2], beta_1[3] +  2*Z[, 2],  beta_1[4]  -  2*Z[, 3],beta_1[5]-2*Z[,4])+.5*rnorm(N)
# coeffs2 <- cbind(beta_2[1]+2*Z[,1], beta_2[2], beta_2[3] +  2*Z[, 2],  beta_2[4] -  2*Z[, 3],beta_2[5]-2*Z[,4])+.5*rnorm(N)
# coeffs3 <- cbind(beta_3[1]+2*Z[,1], beta_3[2], beta_3[3] +  2*Z[, 2],   -  2*Z[, 3],beta_3[5]-2*Z[,4])+.5*rnorm(N)
# coeffs4 <- cbind(beta_4[6]+2*Z[,1], beta_4[7], beta_4[8] +  2*Z[, 2],  beta_4[9]  -  2*Z[, 3],beta_4[10]-2*Z[,4])+.5*rnorm(N)
# coeffs5 <- cbind(beta_5[6]+2*Z[,1], beta_5[7], beta_5[8] +  2*Z[, 2],  beta_5[9]  -  2*Z[, 3],beta_5[10]-2*Z[,4])+.5*rnorm(N)
# coeffs6 <- cbind(beta_6[1]+2*Z[,1], beta_6[2], beta_6[3] +  2*Z[, 2],  beta_6[4]  -  2*Z[, 3],beta_6[5]-2*Z[,4])+.5*rnorm(N)


#
# vProb = cbind((diag(X[, 1:5]%*%t(coeffs1))), (diag(X[, 1:5]%*%t(coeffs2))), (diag(X[, 1:5]%*%t(coeffs3))),(diag(X[, 6:10]%*%t(coeffs4))),(diag(X[, 6:10]%*%t(coeffs5))),(diag(X[, 1:5]%*%t(coeffs6))) )
# y_train=vProb
# #mChoices = t(apply(vProb, 1, FUN=rmultinom, n = 1, size = 1))
# #y_train = apply(mChoices, 1, function(x) which(x==1))
# num_eff=5

theta<-array(0,c(p,K,6))
theta[1,1,1]<-2;theta[3,2,1]<-2;theta[4,3,1]<- -2;theta[5,4,1]<- -2;
theta[1,1,2]<-2;theta[3,2,2]<-2;theta[4,3,2]<- -2;theta[5,4,2]<- -2;
theta[6,1,3]<-2;theta[8,2,3]<-2;theta[9,3,3]<- -2;theta[10,4,3]<- -2;
theta[6,1,4]<-2;theta[8,2,4]<-2;theta[9,3,4]<- -2;theta[10,4,4]<- -2;
theta[11,1,5]<-2;theta[13,2,5]<-2;theta[14,3,5]<- -2;theta[15,4,5]<- -2;
theta[11,1,6]<-2;theta[13,2,6]<-2;theta[14,3,6]<- -2;theta[15,4,6]<- -2

pliable = matrix(0,N,6)
for (e in 1:6) {
  pliable[,e]<-	compute_pliable(X, Z, theta[,,e])

}

esd<-diag(6)
e<-mvrnorm(N,mu=rep(0,6),Sigma=esd)

y_train<-X%*%Beta+pliable+e

y=y_train
colnames(y)<-c(1:6)
colnames(y)<- c( paste("y",1:(ncol(y)),sep = "") )
#num_eff<-4
#coeffs <- cbind(beta[1]+5*Z[,1], beta[2], beta[3] +  3*Z[, 2],  beta[4] *(e -  2*Z[, 3]),beta[5]*(e-2*Z[,4]))



TT=tree.parms(y)
C<-TT$Tree
CW<-TT$Tw


N=nrow(X)
#print(c(N,length(y)))

p=ncol(X)
K=ncol(Z)
D=dim(y)[2]
lambda=rep(0.5,6);alpha=0.5;e.abs=1E-4;e.rel=1E-2;alpha=.2;tol=1E-3;alph=1;rho=5;gg<-c(0.02,0.02);max_it=5000


my_W_hat<-generate.my.w(X=X,Z=Z, quad = TRUE)

svd.w<- svd(my_W_hat)
svd.w$tu<- t(svd.w$u)
svd.w$tv<- t(svd.w$v)



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

  DD1<-rho*(new_I[rr]+1)

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



r_current = y#-model(beta0,theta0,beta=beta_hat, theta, X=W_hat, Z)
b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
beta0<-b$beta0
theta0<-b$theta0

new_y<- y-( matrix(1,N)%*%beta0+Z%*%((theta0)))

XtY <- crossprod((my_W_hat), (new_y))


my_values<-	admm.MADMMplasso(beta0=beta0,theta0=theta0,beta=beta,beta_hat=beta_hat,theta=theta,rho,X,Z,max_it,W_hat=my_W_hat,XtY,y,N,p,K,e.abs, e.rel,alpha, lambda=lambda,alph,svd.w=svd.w,tree = TT,my_print=F,invmat=invmat,V=V,Q=Q,E=E,EE=EE,O=O,P=P,H=H,HH=HH,cv=cv,gg=gg)

beta=my_values$beta;theta=my_values$theta;converge=my_values$converge;beta0=my_values$beta0;theta0=my_values$theta0### iteration
V=my_values$V;Q=my_values$Q;O=my_values$O;P=my_values$P;E=my_values$E;H=my_values$H;beta_hat=my_values$beta_hat; y_hat<-my_values$y_hat
