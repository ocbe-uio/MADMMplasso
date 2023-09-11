# MADMMplasso

Multi variate multi-response 'ADMM' with interaction effects combines the usual squared error loss for the mult-response problem with some penalty terms  to encourage responses that correlate to form groups and also allow for modeling main and interaction effects that exit within the covariates. 

The method can be powperful in situations where one assumes that;
1.	certain factors influence the main covariate seperatly and aims to include these fatcors as modifying varibles to the main covariate. 
2.	There exists some form of grouping within the responses and want to include this information. We assume that the responses form overlapping groups that follows a certain hierarchy. 
A typical example is when one wants to model drug response for multiple drugs and assumes that some of the drugs share certain properties in common, for example drug target and chemical compounds and aims to include this information to improve prediction and also aim to predict which drug could be suitable which patient (disease). The various diseases under study could be the modifying variable. 








Author: Theophilus Asenso, Manuela Zucknick


## Usage

devtools::install_github("ocbe-uio/MADMMplasso")


set.seed(1235)

N = 100 ; p =50;nz=4; K=nz

X <- matrix(rnorm(n = N * p), nrow = N, ncol = p)

mx=colMeans(X)

sx=sqrt(apply(X,2,var))

X=scale(X,mx,sx)

X=matrix(as.numeric(X),N,p)

Z =matrix(rnorm(N*nz),N,nz)

mz=colMeans(Z)

sz=sqrt(apply(Z,2,var))

Z=scale(Z,mz,sz)

beta_1 <- rep(x = 0, times = p)

beta_2<-rep(x = 0, times = p)

beta_3<-rep(x = 0, times = p)

beta_4<-rep(x = 0, times = p)

beta_5<-rep(x = 0, times = p)

beta_6<-rep(x = 0, times = p)




beta_1[1:5] <- c(2, 2, 2, 2,2)

beta_2[1:5]<-c(2, 2, 2, 2,2)

beta_3[6:10]<-c(2, 2, 2, -2,-2)

beta_4[6:10] <- c(2, 2, 2, -2,-2)

beta_5[11:15] <- c(-2,  -2,-2, -2,-2)

beta_6[11:15] <- c(-2, -2, -2, -2,-2)

Beta<-cbind(beta_1,beta_2,beta_3,beta_4,beta_5,beta_6)

colnames(Beta)<-c(1:6)


theta<-array(0,c(p,K,6))

theta[1,1,1]<-2;theta[3,2,1]<-2;theta[4,3,1]<- -2;theta[5,4,1]<- -2;

theta[1,1,2]<-2;theta[3,2,2]<-2;theta[4,3,2]<- -2;theta[5,4,2]<- -2;

theta[6,1,3]<-2;theta[8,2,3]<-2;theta[9,3,3]<- -2;theta[10,4,3]<- -2;

theta[6,1,4]<-2;theta[8,2,4]<-2;theta[9,3,4]<- -2;theta[10,4,4]<- -2;

theta[11,1,5]<-2;theta[13,2,5]<-2;theta[14,3,5]<- -2;theta[15,4,5]<- -2;

theta[11,1,6]<-2;theta[13,2,6]<-2;theta[14,3,6]<- -2;theta[15,4,6]<- -2

library(MASS)

pliable = matrix(0,N,6)
for (e in 1:6) {
  pliable[,e]<-	compute_pliable(X, Z, theta[,,e])
}

esd<-diag(6)

e<-MASS::mvrnorm(N,mu=rep(0,6),Sigma=esd)

y_train<-X%*%Beta+pliable+e

y=y_train

#colnames(y)<-c(1:6)

colnames(y)<- c( paste("y",1:(ncol(y)),sep = "") )

TT=tree.parms(y)

plot(TT$h_clust)



![githubb](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/1a843b46-7154-405c-8db6-cec5b7a0982d)


gg1=matrix(0,2,2)


gg1[1,]<-c(0.02,0.02)

gg1[2,]<-c(0.02,0.02)

nlambda = 50

e.abs=1E-4

e.rel=1E-2

alpha=.5

tol=1E-3


fit <- MADMMplasso(
  X, Z, y, alpha=alpha, my_lambda=NULL,
  lambda_min=0.001, max_it=5000, e.abs=e.abs, e.rel=e.rel, maxgrid=nlambda,
  nlambda=nlambda, rho=5, tree=TT, my_print=FALSE, alph=1, parallel=FALSE,
  pal=1, gg=gg1, tol=tol
)


plot(fit)



![1](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/dbf999bb-d07a-47ff-9a6d-b9ab8cc81dcd)



![2](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/b28f4d07-b634-40be-b303-0335304e1f37)


![3](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/c56e02a8-bff6-4bd7-9e51-6d8fabe5b384)


![4](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/bbdc715a-33bf-4ae9-84c2-cb1852deb860)

![5](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/ff78e14a-00af-4afb-82a2-8b538349b04a)
![6](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/81b87875-511c-4239-ac9d-37c2218688b3)






gg1=fit$gg
cv_admp <- cv.MADMMplasso(
  fit, nfolds=5, X, Z, y, alpha=alpha, lambda=fit$Lambdas, max_it=5000,
  e.abs=e.abs, e.rel=e.rel, nlambda, rho=5, my_print=FALSE, alph=1,
  foldid=NULL, parallel=FALSE, pal=1, gg=gg1, TT=TT, tol=tol
)


plot(cv_admp)






s_ad=which(cv_admp$lambda[,1]==cv_admp$lambda.min)





![cv](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/d5beae83-bf02-475c-999c-d156bd7b9936)





fit$beta[[s_ad]]


![Screenshot 2023-09-11 at 16 03 48](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/0420164e-bfb9-476a-a7ef-3e97bde377e2)







