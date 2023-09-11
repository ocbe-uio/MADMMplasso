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



fit$beta[[s_ad]]


50 x 6 sparse Matrix of class "dgCMatrix"
                                                                                   
 [1,]  2.078508938  2.023153455  0.050872550  0.012623014 -0.063918655  0.022864607
 [2,]  2.061628339  1.985071541  .           -0.085676699 -0.037706542 -0.056261270
 [3,]  1.993691536  2.050836354  .           -0.001136581  .           -0.034209539
 [4,]  1.977787514  2.043220902  0.132675390 -0.059120152  0.147234652  0.066075606
 [5,]  2.161707799  2.045040453  .            .            0.023063060  0.072181837
 [6,]  .            .            1.773454714  1.679624498  .           -0.065286881
 [7,]  .            .            1.755549506  1.850150568  .            0.040566716
 [8,]  .            .            1.974443587  1.764094038  0.009185381  .          
 [9,]  0.020120341 -0.093272501 -1.790483641 -1.873190655  .            .          
[10,]  .            .           -1.951503402 -2.018807953 -0.028990502 -0.056683403
[11,]  .            .            0.166001886  .           -1.756434921 -1.933752819
[12,]  0.105832511  0.035702418  .            .           -2.008163038 -1.809373205
[13,] -0.024291884  .            0.003739564  0.015782320 -1.900751121 -1.684687810
[14,]  0.001481748 -0.037222015  0.013755131  .           -1.751355847 -1.724598326
[15,]  0.005179274  .           -0.064830625 -0.239912483 -1.704251475 -1.756799581
[16,] -0.048395792 -0.056264823  .           -0.047623179 -0.047835878  .          
[17,]  .            .            0.034709139  .            .            .          
[18,]  .            .           -0.005086942  .            .           -0.067470814
[19,]  0.013122300  .            .            .            0.067389844  0.049615630
[20,]  0.009555325  .           -0.312280182  0.018750703  .           -0.012150706
[21,]  .            .            .            .            .            .          
[22,] -0.040967605  0.032717800  0.057185897  0.015493705  .            0.027750422
[23,] -0.018788150  0.059433460  .            .            0.003571943 -0.020855247
[24,]  .            0.050740317  .            .            .            .          
[25,]  .            .            0.028591305  .            .            .          
[26,]  0.034018425  .           -0.108660699  .            .           -0.022508536
[27,] -0.017971226  .            0.166225821 -0.051545535 -0.118122046  .          
[28,] -0.007111410  0.038815800  0.004806171  0.002208667  0.012162400  0.006663189
[29,] -0.025141661  0.175944966 -0.051409257 -0.001184099  0.008716499  .          
[30,]  .           -0.004050746  .           -0.044450455 -0.336872115 -0.089699090
[31,]  0.021612154  .            .            .            .            .          
[32,]  .            .            .            0.024163311  .           -0.037356125
[33,]  0.051773340  .            .            .            .            .          
[34,]  .            .            .            .           -0.045178950  .          
[35,]  .            .            .            .           -0.063895280  .          
[36,]  0.116382895 -0.002486744 -0.002541896  0.088167952 -0.002985550 -0.050974651
[37,]  .            .            .            .            0.018427941  0.053732217
[38,]  .            0.121734040  .            .            .            .          
[39,]  .            .            .            .            .            .          
[40,]  .           -0.001636976  .           -0.080742207 -0.086996943  .          
[41,]  0.045897326  0.080542499  .           -0.040992018  .            .          
[42,]  .           -0.017790812 -0.079020488  .            .            .          
[43,]  .            .            .            .            0.321018378  0.283706761
[44,]  .            .            0.086610844  .            .           -0.144072922
[45,]  0.107412821  .           -0.034959626  .           -0.074577956 -0.091787496
[46,]  0.134476793  .            .            .            .            0.024791837
[47,]  .            .            0.005713297  0.004421865  .            .          
[48,]  0.019199106  0.165041517  0.033969300  0.149092956  .            .          
[49,]  .            .            .            .           -0.003119204  0.117227085
[50,]  .            .           -0.144682422 -0.066361421 -0.004567474  .          
