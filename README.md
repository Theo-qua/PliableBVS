# PliableBVS
A flexible Bayesian variable selection method for modeling interactions


## Usage
```r
devtools::install_github("Theo-qua/PliableBVS")
library(PliableBVS)
set.seed(1)
train_frac = 0.75
val_frac = 0
test_frac = 0.25
N = 5000 ; p =20;nz=4; K=nz
  X <- matrix(rnorm(n = N * p), nrow = N, ncol = p)
  
  mx=colMeans(X)
  
  sx=sqrt(apply(X,2,var))
  X=scale(X,T,F)
  #X=matrix(as.numeric(X),N,p)
  
  #Z <- matrix(rbern(n = N * K,  prob = 0.5), nrow = N, ncol = K)
  Z =matrix(rnorm(N*nz),N,nz)
  mz=colMeans(Z)
  sz=sqrt(apply(Z,2,var))
  
  Z=scale(Z,T,F)
  e=matrix(1,N)
  #X <- matrix(rnorm(n = N * p, mean = 0, sd = 1), nrow = N, ncol = p)
  #Z <- matrix(rbinom(n = N * K, size = 1, prob = 0.5), nrow = N, ncol = K)
  
  e=matrix(1,N)
  
  beta <- Matrix(0,  p,sparse=T)
  beta[1:4] <- c(2, -2, 2, 2)
  coeffs <- cbind(beta[1], beta[2], beta[3] + 2 * Z[, 1], beta[4] * (e - 2 * Z[, 2]))
  theta<-Matrix(0,p,K,sparse = T)
  theta[3,1]<-2;theta[4,2]<- -4#;theta[5,3]<- 2
  #

 
  

  
  
  
 # pliable1<-	compute_pliable(X, Z, theta)


  ##############################################
  num_eff<-4
  #coeffs <- cbind(beta[1]+5*Z[,1], beta[2], beta[3] +  3*Z[, 2],  beta[4] *(e -  2*Z[, 3]),beta[5]*(e-2*Z[,4]))
  signal_to_noise_ratio = 6
 
  mu <-  diag(X[, 1:4] %*% t(coeffs))
  mu <-matrix(mu,N,1)
  
  
  noise = rnorm(N, mean = 0,sd=1)
  
  kk=4
  y_train <- mu + kk*noise
  
  y=matrix(y_train,N,1)
  
  
  
  
  snr = var(as.numeric(mu )) / var(y-as.numeric(mu ))
  
  cat("", fill=T)
  cat(c("snr =",snr),fill=T)
  cat("",fill=T)
  
 
  
  #split training, validation and test sets
 
  smp_size_train = floor(train_frac * nrow(X)) 
  smp_size_val = floor(val_frac * nrow(X))
  train_ind = sort(sample(seq_len(nrow(X)), size = smp_size_train))
  ind_no_train = setdiff(seq_len(nrow(X)), train_ind)
  val_ind = sort(sample(ind_no_train, size = smp_size_val))
  test_ind = setdiff(ind_no_train, val_ind)
  
  colnames(X) = seq(ncol(X))
  colnames(Z) = seq(ncol(Z))
 
  train_x_raw <- X[train_ind, ]
  val_x_raw <- X[val_ind,]
  test_x_raw <- X[test_ind, ]
  
  train_z_raw <- Z[train_ind, ]
  val_z_raw <- Z[val_ind,]
  test_z_raw <- Z[test_ind, ]
  
  
  y_train <- y[train_ind, ]
  val_y <- y[val_ind, ]
  y_test <- y[test_ind, ]
  N_test=length(y_test)
  
  preprocess_values_train_X = preProcess(train_x_raw, method = c("center", "scale"))
  X_train = predict(preprocess_values_train_X, train_x_raw)
  X_test = predict(preprocess_values_train_X, test_x_raw)
  
  preprocess_values_train_Z = preProcess(train_z_raw, method = c("center", "scale"))
  Z_train = predict(preprocess_values_train_Z, train_z_raw)
  Z_test = predict(preprocess_values_train_Z, test_z_raw)
  
 
  
  
  mu_train = mu[train_ind, ]
  mu_test= mu[test_ind, ]
  
  
  err_null <- sum((y_test-mean(y_train))^2)/length(y_test)
  print(err_null)
  
  X=X_train;Z=Z_train;y=y_train
  alpha=0.5
  
  
  system.time(ff<-PliableBVS(Y=y, X,Z,alpha=alpha,family = "gaussian", niter = 5000, burnin = 2000, a_rho = 2, b_rho=3,a_zeta = 2, b_zeta=3,num_update = 50, niter.update =200,burnin.update=100, verbose1 = T,verbose2 = T, lam1=1e-1,lam2=1e-1, rho_prior=TRUE, rho=0.5,zeta=0.5,c2=10^2,v2=10^2, update_tau=TRUE,option.weight.group=FALSE,option.update="global",lambda2_update=NULL) )


plot(ff,type="val",coef_val=c(3))
plot(ff,type="val",coef_val=c(3,1))
plot(ff,type="likelihood",coef_val=c(3,1))
plot(ff,type="ms",coef_val=c(3,1))
plot(ff,type="cont",coef_val=c(3,1))
plot(ff,type="dist",coef_val=c(3,1))
```
![Rplot1](https://github.com/user-attachments/assets/7c6bfce3-d06b-489e-8850-2a16683aeeba)
![Rplot2](https://github.com/user-attachments/assets/52edb2db-abdd-46bc-bc0b-bef761c2b67d)
![Rplot3](https://github.com/user-attachments/assets/18fd2987-f2be-4ff1-ae8b-724c367dcb52)
![Rplot4](https://github.com/user-attachments/assets/b9fe2fd2-6053-45ae-929c-bf9ac3701a45)
![Rplot](https://github.com/user-attachments/assets/8b1a4317-f183-4cb4-98fb-56809d5fe413)
![Rplot6](https://github.com/user-attachments/assets/ef89c878-e66b-4ff4-ac26-9952f9a7b371)


