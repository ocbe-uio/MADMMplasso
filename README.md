# MADMMplasso

Multi variate multi-response 'ADMM' with interaction effects combines the usual squared error loss for the mult-response problem with some penalty terms  to encourage responses that correlate to form groups and also allow for modeling main and interaction effects that exit within the covariates.

The method can be powerful in situations where one assumes that;
1.	certain factors influence the main covariate separately and aims to include these factors as modifying variables to the main covariate.
2.	There exists some form of grouping within the responses and want to include this information. We assume that the responses form overlapping groups that follows a certain hierarchy.
A typical example is when one wants to model drug response for multiple drugs and assumes that some of the drugs share certain properties in common, for example drug target and chemical compounds and aims to include this information to improve prediction and also aim to predict which drug could be suitable for which patient (given a particular disease). The various diseases under study could be the modifying variable.

Author: Theophilus Asenso, Manuela Zucknick


## Usage

```r
devtools::install_github("ocbe-uio/MADMMplasso")
set.seed(1235)
N <- 100; p <- 50; nz <- 4; K <- nz
X <- matrix(rnorm(n = N * p), nrow = N, ncol = p)
mx <- colMeans(X)
sx <- sqrt(apply(X,2,var))
X <- scale(X,mx,sx)
X <- matrix(as.numeric(X),N,p)
Z =matrix(rnorm(N*nz),N,nz)
mz <- colMeans(Z)
sz <- sqrt(apply(Z,2,var))
Z <- scale(Z,mz,sz)
beta_1 <- rep(x = 0, times = p)
beta_2 <- rep(x = 0, times = p)
beta_3 <- rep(x = 0, times = p)
beta_4 <- rep(x = 0, times = p)
beta_5 <- rep(x = 0, times = p)
beta_6 <- rep(x = 0, times = p)

beta_1[1:5] <- c(2, 2, 2, 2,2)
beta_2[1:5] <- c(2, 2, 2, 2,2)
beta_3[6:10] <- c(2, 2, 2, -2,-2)
beta_4[6:10] <- c(2, 2, 2, -2,-2)
beta_5[11:15] <- c(-2,  -2,-2, -2,-2)
beta_6[11:15] <- c(-2, -2, -2, -2,-2)

Beta<-cbind(beta_1,beta_2,beta_3,beta_4,beta_5,beta_6)

colnames(Beta) <- c(1:6)

theta <- array(0,c(p,K,6))
theta[1,1,1] <- 2; theta[3,2,1] <- 2; theta[4,3,1] <- -2; theta[5,4,1] <- -2;
theta[1,1,2] <- 2; theta[3,2,2] <- 2; theta[4,3,2] <- -2; theta[5,4,2] <- -2;
theta[6,1,3] <- 2; theta[8,2,3] <- 2; theta[9,3,3] <- -2; theta[10,4,3] <- -2;
theta[6,1,4] <- 2; theta[8,2,4] <- 2; theta[9,3,4] <- -2; theta[10,4,4] <- -2;
theta[11,1,5] <- 2; theta[13,2,5] <- 2; theta[14,3,5] <- -2; theta[15,4,5] <- -2;
theta[11,1,6] <- 2; theta[13,2,6] <- 2; theta[14,3,6] <- -2; theta[15,4,6] <- -2

library(MASS)
pliable <- matrix(0,N,6)
for (e in 1:6) {
  pliable[,e] <- compute_pliable(X, Z, theta[,,e])
}
esd <- diag(6)
e <- MASS::mvrnorm(N,mu=rep(0,6),Sigma=esd)
y_train <- X %*% Beta + pliable + e
y <- y_train
colnames(y) <- c( paste("y",1:(ncol(y)),sep = "") )
TT <- tree_parms(y)
plot(TT$h_clust)
```

![github](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/1a843b46-7154-405c-8db6-cec5b7a0982d)

```r
gg1 <- matrix(0,2,2)
gg1[1,] <- c(0.02,0.02)
gg1[2,] <- c(0.2,0.2)
nlambda <- 50
e.abs <- 1E-4
e.rel <- 1E-2
alpha <- .5
tol <- 1E-3

fit <- MADMMplasso(
  X, Z, y, alpha=alpha, my_lambda=NULL,
  lambda_min=0.001, max_it=5000, e.abs=e.abs, e.rel=e.rel, maxgrid=nlambda,
  nlambda=nlambda, rho=5, tree=TT, my_print=FALSE, alph=1, gg=gg1, tol=tol
)

plot(fit)
```

![1](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/b8841ba1-aac6-4539-9924-70c70accddd9)
![2](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/c2e4bfcf-22c8-49a7-bf99-07ddb436437b)
![3](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/b319ad79-71bf-4de2-9d9e-457f50393a1e)
![4](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/34d8d6e1-c912-4654-a497-4bade67d5ee1)
![5](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/fe375fff-51e2-4b49-9520-f7cbcaec6bbb)
![6](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/c4c46d9b-3cd3-4c55-95d1-abbb59405422)

```r
gg1 <- fit$gg

cv_admp <- cv_MADMMplasso(
  fit, nfolds=5, X, Z, y, alpha=alpha, lambda=fit$Lambdas, max_it=5000,
  e.abs=e.abs, e.rel=e.rel, nlambda, rho=5, my_print=FALSE, alph=1,
  foldid=NULL, gg=gg1, TT=TT, tol=tol
)

plot(cv_admp)
```

![cv](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/0118f157-dd7a-4387-88f9-f0e18434d59d)

```r
s_ad <- which(cv_admp$lambda[,1]==cv_admp$lambda.min)
fit$beta[[s_ad]]
```

![Screenshot 2023-09-11 at 16 25 59](https://github.com/ocbe-uio/MADMMplasso/assets/85598983/f762b9e1-9212-43c7-a21c-b83a9a48662f)
