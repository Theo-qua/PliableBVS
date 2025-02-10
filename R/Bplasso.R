#' @title Compute the interaction part of the model.
#' @description  Compute the interaction part of the model.
#' @param X N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical variables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param theta    theta coefficients for a single response ncol(X) by ncol(Z)
#' @return a vector of length N of the calculated interaction term for a single response
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


#
# library(pracma)
#
# reg<-function(r,Z){
#   K=ncol(Z)
#   N=nrow(Z)
#   my_one<-matrix(1,nrow(Z))
#   my_w=data.frame(Z)
#   my_w<-as.matrix(my_w)
#   my_inv<-pinv( (t(my_w)%*%my_w)/N )
#   #my_res<-my_inv%*%( (t(my_w)%*%r)/N )
#   # new<- lm(r~1,na.action=na.exclude)
#   #beta0<-matrix(my_res[(K+1)])
#   new1<- lm(r~Z,singular.ok = TRUE)
#   beta0<-matrix(new1$coefficients[1])
#
#   theta0<- matrix(new1$coefficients[-1])
#   #theta0<- matrix(my_res[c(1:(K))])
#   return(list(beta0,theta0))
# }


#' @title Compute the entire regression response model.
#' @description  Compute the entire regression response model.
#' @param X N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical variables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param theta    theta coefficients for a single response ncol(X) by ncol(Z)
#' @param beta a vectore of estimated beta coefficients each having a length   ncol(X)
#' @param theta0  (main effect) coefficients for Z ncol(Z)
#' @param beta0   beta_0 coefficients
#' @return a vector of length N of the  response


#' @export
model<-function(beta0, theta0, beta, theta, X, Z){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)


  intercepts = as.numeric(beta0)+Z%*%(matrix(theta0))
  shared_model = X%*%matrix(beta)
  pliable = compute_pliable(X, Z, theta)
  return(intercepts+  shared_model +pliable)
}


errfun.binomial=function(y,yhat,w=rep(1,length(y))){
  prob_min = 1e-05
  prob_max = 1 - prob_min
  predmat = pmin(pmax(yhat, prob_min), prob_max)

  -2*w*(y*log(predmat)+(1-y)*log(1-predmat))
}



#' @title Fit a a Bayesian pliable lasso model over some MCMC iterations
#' @description Fit a a Bayesian pliable lasso model over some MCMC iterations
#' @param X  N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical variables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param y N by D matrix  of responses. The X and Z variables are centered in the function. We recommend that X and Z also be standardized before the call
#' @param family response type- either "gaussian", "binomial". In the binomial case, y should be 0s and 1s.
#' @return a list of all the posterior objects
#' @export
Bpliable = function(Y, X,Z,alpha=0.5,family = c("gaussian", "binomial"), niter = 10000, burnin = 5000, a_rho=1, b_rho=1,a_zeta=1, b_zeta=1,num_update = 100, niter.update =100,burnin.update=50, verbose1 = FALSE,verbose2 = FALSE, lam1=1e-1,lam2=1e-1, rho_prior=TRUE, rho=0.5,zeta=0.5,c2=10^2,v2=1e-1, update_tau=TRUE,option.weight.group=FALSE,option.update="global",lambda2_update=NULL,nethod){
  this.call = match.call()
  if(family=="gaussian"){
    cat("using family=",family,",so the loss function is Gaussian")
    fit=Bpliable_gs(Y, X,Z,alpha=alpha, niter = niter, burnin = burnin, a_rho=a_rho, b_rho=b_rho,a_zeta=a_zeta, b_zeta=b_zeta,num_update = num_update, niter.update =niter.update,burnin.update=burnin.update, verbose1 = verbose1,verbose2 = verbose2, lam1=lam1,lam2=lam2, rho_prior=rho_prior, rho=rho,zeta=zeta,c2=c2,v2=v2, update_tau=update_tau,option.update=option.update,lambda2_update=lambda2_update)
  }else{
    cat("using family=",family,",so the loss function is logistic")
    fit=Bpliable_lr(Y, X,Z,alpha=alpha, niter = niter, burnin = burnin, a_rho=a_rho, b_rho=b_rho,a_zeta=a_zeta, b_zeta=b_zeta,num_update = num_update, niter.update =niter.update,burnin.update=burnin.update, verbose1 = verbose1,verbose2 = verbose2, lam1=lam1,lam2=lam2, rho_prior=rho_prior, rho=rho,zeta=zeta,c2=c2,v2=v2, update_tau=update_tau,option.update=option.update,lambda2_update=lambda2_update)

  }

  fit$call=this.call
  class(fit) <- "Bpliable"

  return(fit)
}


Bpliable_gs = function(Y, X,Z,alpha=0.5, niter = 10000, burnin = 5000, a_rho=1, b_rho=1,a_zeta=1, b_zeta=1,num_update = 100, niter.update =100,burnin.update=50, verbose1 = FALSE,verbose2 = FALSE, lam1=1e-1,lam2=1e-1, rho_prior=TRUE, rho=0.5,zeta=0.5,c2=10^2,v2=1e-1, update_tau=TRUE,option.weight.group=FALSE,option.update="global",lambda2_update=NULL)
{
  ####################################
  # Create and Initialize parameters #
  ####################################
  n = length(Y)
  p = dim(X)[2]
  K=dim(Z)[2]

  # Initialize parameters
  tau_beta2 = rep(1, p)
  tau_theta2 = matrix(1, p,K)
  tau_theta0=1
  sigma2 = 1
  eta = rep(0, p)
  gamma=matrix(0,p,K  )

  beta = rep(1, p)
  theta=matrix(1,p,K  )


  beta0<-1
   theta0<-rep(1,K)


  Q = rep(0,p)
  R<-matrix(0,p,K)


  #my_W_hat<-generate.my.w(X=X,Z=Z, quad = TRUE)
  ##################################
  # Compute lambda2 via EM         #
  ##################################
  if (update_tau==TRUE){
    fit_for_lambda2 = Bpliable_EM_lambda_gs(Y, X,Z,alpha, num_update = num_update, niter = niter.update,burnin =burnin.update , a_rho=a_rho, b_rho=b_rho,a_zeta=a_zeta, b_zeta=b_zeta,lam1=lam1,lam2=lam2,c2=c2,v2=v2,option.update=option.update,option.weight.group=option.weight.group,verbose = verbose2)
    lambda2 = apply(fit_for_lambda2$lambda2_path,2,tail,1)}else
    {
      lambda2 <- rep(lambda2_update,p)
    }

#print(lambda2)
  #lambda2 = rep(.01,p)

  #print(c("lambda2",option.update))


mean_Y<-mean(Y)
Y_center<-Y-mean(Y)

  #####################
  # The Gibbs Sampler #
  #####################

  # burnin
  coef = array(0, dim=c(niter-burnin,p ))
coef_theta =array(0, dim=c(niter-burnin,p,K ))

coef_all = array(0, dim=c(niter, p))
coef_theta_all =array(0, dim=c(niter,p,K))

Q_all = array(0, dim=c(niter, p))
R_all =array(0, dim=c(niter,p,K))

coef_theta1 =list()
  coef_tau = array(0, dim=c(niter, p))

  coef_beta0 = array(0, dim=c(niter-burnin, 1))
  coef_theta0 = array(0, dim=c(niter-burnin, K))

 likelihood = array(0, dim=c(niter, 1))
  count_coeff=1
  for (iter in 1:niter)
  {
    # print the current number of iteration
    if ( isTRUE(verbose1 == TRUE & iter %% 10==0)==TRUE | isTRUE(verbose1 == TRUE & iter ==1)==TRUE) {print(iter)}

   #  r_current = Y-model(beta0, theta0, beta=beta, theta, X=X, Z)
   # bb = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
   #  beta0<-as.numeric(unlist(bb[1]))
   #  theta0<-matrix(unlist(bb[2]))
   #
   #

   #
    Y_bar<-Y-as.numeric(beta0)-Z%*%(matrix(theta0))-X%*%beta-compute_pliable(X, Z, theta)

    #Y_bar<-Y_center-Z%*%(matrix(theta0))-X%*%beta-compute_pliable(X, Z, theta)

    Y_bar1=Y_bar+as.numeric(beta0)
    beta01=beta0

    beta0<-rmnorm(1,mean=sum(Y_bar1)/(n*(1+sigma2/c2)),vcov=1/(n*(1+sigma2/c2)))

    Y_bar<-Y_bar-as.numeric(beta0)+as.numeric(beta01)

    theta01=theta0
    Y_bar1=Y_bar+Z%*%(matrix(theta0))
    #f_theta0<-(t(Z)%*%Z)
   #  f_theta0<-(1/tau_theta0+t(Z)%*%Z)
   #  f_theta0_inverse<-solve(f_theta0)
   # YtZ<-t(Y_bar1)%*%Z
   # mean_theta0<-f_theta0_inverse%*%t(YtZ)
   #
   #  theta0<-rmnorm(1,mean=mean_theta0,vcov=sigma2*f_theta0_inverse)
   #

    f_theta0<-(sigma2/v2+t(Z)%*%Z)
    f_theta0_inverse<-solve(f_theta0)
    YtZ<-t(Y_bar1)%*%Z
    mean_theta0<-f_theta0_inverse%*%t(YtZ)

    theta0<-rmnorm(1,mean=mean_theta0,vcov=f_theta0_inverse)





    Y_bar<- Y_bar-Z%*%(matrix(theta0))+Z%*%(matrix(theta01))

    intercepts=as.numeric(beta0)+Z%*%(matrix(theta0))

    # Update beta's
    for(j in 1:p)
    {
      beta_j=beta[j]

      #Y_bar<-Y- intercepts-X%*%beta-compute_pliable(X, Z, theta)

      f1 = t(X[,j])%*%(Y_bar+X[,j]*beta[j])
      f2 = t(X[,j])%*%X[,j]+1/tau_beta2[j]
      f2_inverse = 1/f2#solve(f2)
      mu = f2_inverse %*% f1
      ### Main part
      eta[j] = rho/(rho+(1-rho)*(tau_beta2[j])^(-1/2)*(f2)^(-1/2)*exp(t(f1)%*%mu/(2*sigma2)))
      maxf <- max(f2)
      trythis <- (-1/2)*log(tau_beta2[j]) + (-1/2)*log((f2/maxf)) + (-dim(f2)[1]/2)*log(maxf) + t(f1)%*%mu/(2*sigma2)
      eta[j] = rho/(rho+(1-rho)*exp(trythis))

      if(runif(1)<eta[j])
      {
        beta[j] = 0
        Q[j] = 0

        theta[j,]<-0
        R[j,]<-0


      }
      else
      {
        beta[j] = rmnorm(1, mean=mu, vcov=sigma2*f2_inverse)
        Q[j] = 1
        Y_bar=Y_bar+X[,j]*beta_j-X[,j]*beta[j]
        #print(beta[j])
        ## theta_j update for non zero beta_j
        Y_bar_j<-Y_bar#+X[,j]*beta_j-X[,j]*beta[j]
        for (k in 1:K) {
          theta_jk=theta[j,k]

          f1 = t(X[,j]*Z[,k])%*%(Y_bar_j+X[,j]*Z[,k]*theta[j,k])
          f2 = t(X[,j]*Z[,k])%*%(X[,j]*Z[,k])+1/tau_theta2[j,k]
          f2_inverse = 1/f2#solve(f2)
          mu = f2_inverse %*% f1
          ### Main part
          gamma[j,k] = zeta/(zeta+(1-zeta)*(tau_theta2[j,k])^(-1/2)*(f2)^(-1/2)*exp(t(f1)%*%mu/(2*sigma2)))
          maxf <- max(f2)
          trythis <- (-1/2)*log(tau_theta2[j,k]) + (-1/2)*log((f2/maxf)) + (-dim(f2)[1]/2)*log(maxf) + t(f1)%*%mu/(2*sigma2)
          gamma[j,k] = zeta/(zeta+(1-zeta)*exp(trythis))

          if(runif(1)<gamma[j,k])
          {


            theta[j,k]<-0
            R[j,k]<-0


          }
          else
          {

            theta[j,k] = rmnorm(1, mean=mu, vcov=sigma2*f2_inverse)
            R[j,k] = 1

            Y_bar_j=Y_bar_j+X[,j]*Z[,k]*theta_jk-X[,j]*Z[,k]*theta[j,k]



          }

        }

        Y_bar=Y_bar_j
        #print(c(iter,j,theta[j,]))
      }


    }

    # Update tau2's and tau_theta2's



    if(update_tau) {
      tau_theta0=1/rig(1, mean=sqrt( ( v2*sigma2)/sum(theta0^2)), scale = 1/( v2 ) )


        for(j in 1:p)
        {
          if(Q[j]==0){tau_beta2[j] = rgamma(1, shape=1, rate=( ((1-alpha)^2)*lambda2[j])/2)}
          else{tau_beta2[j] = 1/rig(1, mean=sqrt( ( ((1-alpha)^2)*lambda2[j]*sigma2)/sum(beta[j]^2)), scale = 1/( ((1-alpha)^2)*lambda2[j]))}

          for (k in 1:K) {


          if(R[j,k]==0){tau_theta2[j,k] = rgamma(1, shape=1, rate=( ((alpha)^2)*lambda2[j])/2)}
          else{tau_theta2[j,k] = 1/rig(1, mean=sqrt( ( ((alpha)^2)*lambda2[j]*sigma2)/sum(theta[j,k]^2)), scale = 1/( ((alpha)^2)*lambda2[j]))}

          }

        }

    }


    # Update sigma2
    Y_bar<-Y- intercepts-X%*%beta-compute_pliable(X, Z, theta)
    s=0
    ss=0
    for(i in 1:p)
    {
      s = s + sum(beta[j]^2)/tau_beta2[j]

      for (k in 1:K) {
        ss = ss + sum(theta[j,k]^2)/tau_theta2[j,k]
      }
    }
    beta_vec = c(beta)
    beta0_vec=c(beta0)
    theta0_vec=c(theta0)


    coef_all[iter,] = beta
    Q_all[iter,]=Q

    for (nn in 1:K) {
      coef_theta_all[iter,,nn] =theta[,nn]
      R_all[iter,,nn] =R[,nn]
    }


    if(iter > burnin){
      coef[iter-burnin,] = beta_vec

    for (nn in 1:K) {
      coef_theta[iter-burnin,,nn]=theta[,nn]
    }

    coef_theta1[[count_coeff]]=theta
    count_coeff=count_coeff+1
    coef_beta0[iter-burnin,] = beta0_vec
    coef_theta0[iter-burnin,] = theta0_vec
    }
    sigma2 = rinvgamma(1, shape=(n)/2 + sum(Q)/2+sum(R)/2 + lam1,
                       scale=( t(Y_bar)%*%Y_bar+s+ss)/2 + lam2)

    #sigma2 = rinvgamma(1, shape=(n)/2+(K)/2 + sum(Q)/2+sum(R)/2 + lam1,
     #                  scale=( t(Y_bar)%*%Y_bar+s+ss+sum(theta0^2)/tau_theta0 )/2 + lam2)
    #sigma2 = rinvgamma(1, shape=(n)/2+(K)/2 + sum(Q)/2+sum(R)/2 + lam1,
    #                   scale=( t(Y_bar)%*%Y_bar+s+ss+sum(theta0))/2 + gamma)

  # sigma2 = rinvgamma(1, shape=(n)/2+ sum(Q)/2 +sum(R)/2 + lam1,
  #                    scale=( t(Y_bar)%*%Y_bar+s+ss+beta0+sum(theta0))/2 + gamma)


    # Update pi
    if(rho_prior==TRUE){
      rho = rbeta(1, shape1=a_rho+p-sum(Q), shape2=b_rho+sum(Q))
      zeta=rbeta(1, shape1=a_zeta+p*K-sum(R), shape2=b_zeta+sum(R))

      #rho = rbeta(1, shape1=a+sum(Q), shape2=b+p-sum(Q))
      #zeta=rbeta(1, shape1=a+sum(R), shape2=b+p*K-sum(R))
    }

    likelihood[iter,] = (( ( ( (2*pi*sigma2))^(-n/2) ) ) *exp( (-1/(2*sigma2))*(sum(Y_bar^2))) )
  }

  # output the posterior mean and median as our estimator
  pos_mean_beta = apply(coef, 2, mean)
  pos_median_beta = apply(coef, 2, median)
  pos_mean_theta = (Reduce("+",coef_theta1))/(niter-burnin)
  pos_median_theta<-matrix(0,p,K)
  pos_mean_theta<-matrix(0,p,K)

  for (j in 1:p) {
    for (k in 1:K) {
      pos_median_theta[j,k]= median(coef_theta[,j,k])
      pos_mean_theta[j,k]= mean(coef_theta[,j,k])
    }
  }



  pos_mean_beta0 = apply(coef_beta0, 2, mean)
  pos_median_beta0 = apply(coef_beta0, 2, median)

  pos_mean_theta0 = apply(coef_theta0, 2, mean)
  pos_median_theta0 = apply(coef_theta0, 2, median)

  #num_iter=niter-burnin
  new_beta<-matrix(0,p)
  new_theta<-matrix(0,p,K)
  new_Q=Q_all[-c(1:burnin),]
  new_R=R_all[-c(1:burnin),,]
  for (nn in 1:p) {
    if( (sum(new_Q[,nn])/dim(new_Q)[1])>0.5 ){
      new_beta[nn]=sum(coef[,nn])/sum(new_Q[,nn])
    }

    for (n in 1:K) {

      if( (sum(new_R[,nn,n])/dim(new_R)[1])>0.5 ){
        new_theta[nn,n]=sum(coef_theta[,nn,n])/sum(new_R[,nn,n])
      }
    }
  }


  list(pos_mpm_beta=new_beta,pos_mean_beta = pos_mean_beta, pos_median_beta = pos_median_beta, coef = coef,pos_mean_beta0 = pos_mean_beta0, pos_median_beta0 = pos_median_beta0, coef_beta0 = coef_beta0,pos_mean_theta0 = pos_mean_theta0, pos_median_theta0 = pos_median_theta0, coef_theta0 = coef_theta0,pos_mpm_theta=new_theta ,pos_mean_theta = pos_mean_theta, pos_median_theta = pos_median_theta, coef_theta = coef_theta,coef_beta_all=coef_all,coef_theta_all=coef_theta_all,EM_fit=fit_for_lambda2,lambda2=lambda2,Q=Q_all,R=R_all,Likelihood=likelihood)



}

###########################################################


Bpliable_EM_lambda_gs = function(Y, X,Z,alpha=0.5, num_update = 100, niter = 100,burnin = 50, a_rho=1, b_rho=1,a_zeta=1, b_zeta=1,verbose = FALSE, lam1=1e-1, lam2=1e-1, rho_prior=TRUE, rho=0.5,zeta=0.5,c2=100,v2=1e-1,option.update="global",option.weight.group=FALSE)
{
  ####################################
  # Create and Initialize parameters #
  ####################################
  n = length(Y)
  p = dim(X)[2]
  K=dim(Z)[2]

  # initialize parameters
  tau_beta2 = rep(1, p)
  tau_theta2 = matrix(1, p,K)
  tau_theta0=1

  sigma2 = 1
  lambda2 = 1
  matlambda2 = rep(1,p)
  lambda2_path = rep(-1, num_update)
  matlambda2_path = matrix(-1,ncol=p,nrow=num_update)

  coef_path = array(0, dim=c(p, num_update))
  coef_median_path = array(0, dim=c(p, num_update))
  coef_theta_path =array(0, dim=c(p,K, num_update))
  coef_theta_median_path =array(0, dim=c(p,K, num_update))


  tau2_each_update_path = array(0, dim=c(p, num_update))
  tau2_each_update_median_path = array(0, dim=c(p, num_update))
  tau2_each_update_theta_path = list()

  coef_beta0_path = array(0, dim=c(1, num_update))
  coef_beta0_median_path = array(0, dim=c(1, num_update))
  coef_theta0_path = array(0, dim=c(K, num_update))
  coef_theta0_median_path = array(0, dim=c(K, num_update))



  eta = rep(0, p)
  gamma=matrix(0,p,K  )

  beta = rep(0, p)
  theta=matrix(0,p,K  )
  beta0<-0
  theta0<-rep(0,K)

  Q = rep(0,p)
  R<-matrix(0,p,K)

  mean_Y<-mean(Y)
  Y_center<-Y-mean(Y)
  #####################
  # The Gibbs Sampler #
  #####################

  for (update in 1:num_update) {
    # print(c("updadte=",update))
    coef = array(0, dim=c(p, niter))
    coef_theta =array(0, dim=c(p,K, niter))
    coef_theta1 =list()

    tau2_each_update = array(0, dim=c(p, niter))
    tau2_each_update_theta = list()

    coef_beta0 = array(0, dim=c(1, niter))
    coef_theta0 = array(0, dim=c(K, niter))
    count_coeff=1
    count_tau2=1
    for (iter in 1:niter)    {
      # print the current number of iteration
      if (verbose == TRUE) {print(c(update,iter))}
      # r_current = Y-model(beta0, theta0, beta=beta, theta, X=X, Z)
      # bb = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
      # beta0<-as.numeric(unlist(bb[1]))
      # theta0<-matrix(unlist(bb[2]))
      #


      Y_bar<-Y-as.numeric(beta0)-Z%*%(matrix(theta0))-X%*%beta-compute_pliable(X, Z, theta)

     # Y_bar<-Y_center-Z%*%(matrix(theta0))-X%*%beta-compute_pliable(X, Z, theta)

      Y_bar1=Y_bar+as.numeric(beta0)
      beta01=beta0

      beta0<-rmnorm(1,mean=sum(Y_bar1)/(n*(1+sigma2/c2)),vcov=1/(n*(1+sigma2/c2)))

      Y_bar<-Y_bar-as.numeric(beta0)+as.numeric(beta01)

      theta01=theta0
      Y_bar1=Y_bar+Z%*%(matrix(theta0))
      #f_theta0<-(1+t(Z)%*%Z)
      # f_theta0<-(t(Z)%*%Z)
      # f_theta0_inverse<-solve(f_theta0)
      # YtZ<-t(Y_bar1)%*%Z
      # mean_theta0<-f_theta0_inverse%*%t(YtZ)
      #
      # theta0<-rmnorm(1,mean=mean_theta0,vcov=sigma2*f_theta0_inverse)
      #
      f_theta0<-(sigma2/v2+t(Z)%*%Z)
      f_theta0_inverse<-solve(f_theta0)
      YtZ<-t(Y_bar1)%*%Z
      mean_theta0<-f_theta0_inverse%*%t(YtZ)

      theta0<-rmnorm(1,mean=mean_theta0,vcov=f_theta0_inverse)


      Y_bar<- Y_bar-Z%*%(matrix(theta0))+Z%*%(matrix(theta01))

      intercepts= as.numeric(beta0)+Z%*%(matrix(theta0))

      # Update beta's
      for(j in 1:p)
      {
        beta_j=beta[j]
        #Y_bar<-Y- intercepts-X%*%beta-compute_pliable(X, Z, theta)
        f1 = t(X[,j])%*%(Y_bar+X[,j]*beta[j])
        f2 = t(X[,j])%*%X[,j]+1/tau_beta2[j]
        f2_inverse = 1/f2#solve(f2)
        mu = f2_inverse %*% f1

        eta[j] = rho/(rho+(1-rho)*(tau_beta2[j])^(-1/2)*(f2)^(-1/2)*exp(t(f1)%*%mu/(2*sigma2)))
        maxf <- max(f2)
        trythis <- (-1/2)*log(tau_beta2[j]) + (-1/2)*log((f2/maxf)) + (-dim(f2)[1]/2)*log(maxf) + t(f1)%*%mu/(2*sigma2)
        eta[j] = rho/(rho+(1-rho)*exp(trythis))


        if(runif(1)<eta[j])
        {
          beta[j] = 0
          Q[j] = 0



          theta[j,]<-0
          R[j,]<-0
        }else
        {
          beta[j] = rmnorm(1, mean=mu, vcov=sigma2*f2_inverse)
          Q[j] = 1

          Y_bar=Y_bar+X[,j]*beta_j-X[,j]*beta[j]

         # print( beta[j])
         # Y_bar_j<-Y- intercepts-X%*%beta-compute_pliable(X, Z, theta_j)
          Y_bar_j<-Y_bar#+X[,j]*beta_j-X[,j]*beta[j]
          ## theta_j update for non zero beta_j

          for (k in 1:K) {
            theta_jk=theta[j,k]

            f1 = t(X[,j]*Z[,k])%*%(Y_bar_j+X[,j]*Z[,k]*theta[j,k])
            f2 = t(X[,j]*Z[,k])%*%(X[,j]*Z[,k])+1/tau_theta2[j,k]
            f2_inverse = 1/f2#solve(f2)
            mu = f2_inverse %*% f1
            ### Main part
            gamma[j,k] = zeta/(zeta+(1-zeta)*(tau_theta2[j,k])^(-1/2)*(f2)^(-1/2)*exp(t(f1)%*%mu/(2*sigma2)))
            maxf <- max(f2)
            trythis <- (-1/2)*log(tau_theta2[j,k]) + (-1/2)*log((f2/maxf)) + (-dim(f2)[1]/2)*log(maxf) + t(f1)%*%mu/(2*sigma2)
            gamma[j,k] = zeta/(zeta+(1-zeta)*exp(trythis))

            if(runif(1)<gamma[j,k])
            {


              theta[j,k]<-0
              R[j,k]<-0

            }
            else
            {

              theta[j,k] = rmnorm(1, mean=mu, vcov=sigma2*f2_inverse)
              R[j,k] = 1

              Y_bar_j=Y_bar_j+X[,j]*Z[,k]*theta_jk-X[,j]*Z[,k]*theta[j,k]



            }

          }

          Y_bar=Y_bar_j
        }
      }

      # Update tau2's

      for(j in 1:p)
      {
        if(Q[j]==0){tau_beta2[j] = rgamma(1, shape=1, rate=((1-alpha)^2)*matlambda2[j]/2)}
        else{tau_beta2[j] = 1/rig(1, mean=sqrt(( ((1-alpha)^2)*matlambda2[j]*sigma2)/sum(beta[j]^2)), scale = 1/( ((1-alpha)^2)*matlambda2[j]))}


        for (k in 1:K) {


          if(R[j,k]==0){tau_theta2[j,k] = rgamma(1, shape=1, rate=( ((alpha)^2)*matlambda2[j])/2)}
          else{tau_theta2[j,k] = 1/rig(1, mean=sqrt( ( ((alpha)^2)*matlambda2[j]*sigma2)/sum(theta[j,k]^2)), scale = 1/( ((alpha)^2)*matlambda2[j]))}

        }
      }

      tau_theta0=1/rig(1, mean=sqrt( ( v2*sigma2)/sum(theta0^2)), scale = 1/( v2 ) )



      tau2_each_update[,iter] = tau_beta2
      tau2_each_update_theta[[count_tau2]]=tau_theta2
      count_tau2=count_tau2+1

      # Update sigma2
      Y_bar<-Y-  intercepts-X%*%beta-compute_pliable(X, Z, theta)
      s=0
      ss=0
      for(i in 1:p)
      {
        s = s + sum(beta[j]^2)/tau_beta2[j]
        for (k in 1:K) {
          ss = ss + sum(theta[j,k]^2)/tau_theta2[j,k]
        }
      }
      beta_vec = c(beta)
      beta0_vec=c(beta0)
      theta0_vec=c(theta0)




      if(iter > burnin){
        coef[,iter-burnin] = beta_vec
      coef_theta[,,iter-burnin]=theta
      coef_theta1[[count_coeff]]=theta
      count_coeff=count_coeff+1
      coef_beta0[,iter-burnin] = beta0_vec
      coef_theta0[,iter-burnin] = theta0_vec
      }

      #print(c(s,ss))
     # sigma2 = rinvgamma(1, shape=(n)/2+(K)/2 + sum(Q)/2+sum(R)/2 + lam1,
      #                   scale=( t(Y_bar)%*%Y_bar+s+ss+sum(theta0^2)/tau_theta0 )/2 + lam2)

      sigma2 = rinvgamma(1, shape=(n)/2 + sum(Q)/2+sum(R)/2 + lam1,
                         scale=( t(Y_bar)%*%Y_bar+s+ss)/2 + lam2)


     # sigma2 = rinvgamma(1, shape=(n)/2+(K+1)/2 + sum(Q)/2+sum(R)/2 + lam1,
      #                   scale=( t(Y_bar)%*%Y_bar+s+ss+(beta0^2)/c2+sum(theta0^2)/tau_theta0 )/2 + lam2)

     # sigma2 = rinvgamma(1, shape=(n)/2+(K)/2 + sum(Q)/2+sum(R)/2 + lam1,
      #                   scale=( t(Y_bar)%*%Y_bar+s+ss+sum(theta0))/2 + gamma)

      #sigma2 = rinvgamma(1, shape=(n)/2 + sum(Q)/2+sum(R)/2 + lam1,
     #                    scale=( t(Y_bar)%*%Y_bar+s+ss+beta0+sum(theta0))/2 + gamma)

      # Update phi and zeta
      if(rho_prior==TRUE){
        # rho = rbeta(1, shape1=a+sum(Q), shape2=b+p-sum(Q))
        # zeta=rbeta(1, shape1=a+sum(R), shape2=b+p*K-sum(R))

        rho = rbeta(1, shape1=a_rho+p-sum(Q), shape2=b_rho+sum(Q))
        zeta=rbeta(1, shape1=a_zeta+p*K-sum(R), shape2=b_zeta+sum(R))
      }
    }

    # Update lambda
    tau2_mean = apply(tau2_each_update, 1, mean)
    tau2_mean_theta=(Reduce("+",tau2_each_update_theta))/(niter)

   matlambda2 = (2*p +2*p*K) / ( ( (1-alpha)^2)*sum(tau2_mean) +(alpha^2)*sum(tau2_mean_theta) )

   lambda2 = (2*p +2*p*K) / ( ( (1-alpha)^2)*sum(tau2_mean) +(alpha^2)*sum(tau2_mean_theta))


    #matlambda2 = (2*p ) / ( ( (1-alpha)^2)*sum(tau2_mean)  )

   #lambda2 = (2*p) / ( ( (1-alpha)^2)*sum(tau2_mean) )


    if(option.update=="global") matlambda2 <- rep(lambda2,p)

    matlambda2_path[update,] = matlambda2




    # output the posterior mean and median as our estimator
    coef_path[,update] = apply(coef, 1, mean)
    coef_median_path[,update] = apply(coef, 1, median)
    pos_mean_theta = (Reduce("+",coef_theta1))/(niter)
    pos_median_theta<-matrix(0,p,K)
    pos_mean_theta<-matrix(0,p,K)

    for (j in 1:p) {
      for (k in 1:K) {
        coef_theta_median_path[j,k,update]= median(coef_theta[j,k,])
        coef_theta_path[j,k,update]= mean(coef_theta[j,k,])
      }
    }




    coef_beta0_path[update] = apply(coef_beta0, 1, mean)
    coef_beta0_median_path[update] = apply(coef_beta0, 1, median)

    coef_theta0_path[,update] = apply(coef_theta0, 1, mean)
    coef_theta0_median_path[,update] = apply(coef_theta0, 1, median)


  }



  list(pos_mean_beta = coef_path, pos_median_beta = coef_median_path, pos_mean_beta0 = coef_beta0_path, pos_median_beta0 = coef_beta0_median_path,pos_mean_theta0 = coef_theta0_path, pos_median_theta0 = coef_theta0_median_path,pos_mean_theta = coef_theta_path, pos_median_theta = coef_theta_median_path, lambda2_path = matlambda2_path)
}


#' @export
plot.Bpliable <- function(x, type=c("likelihood","dist","val","cont","ms"),coef_val=c(1,1),...) {
  ff <- x
  if (type=="likelihood") {
    plot(log(ff$Likelihood),type = "l",ylab = "log(likelihood)",xlab = "iterations")
  }else if(type=="ms"){
    model_size<-c()
    for (i in 1:nrow(ff$coef_beta_all)) {
      beta_nz<-sum(ff$coef_beta_all[i,]!=0)

      theta_nz<-sum(ff$coef_theta_all[i,,]!=0)
      model_size=c(model_size,sum(beta_nz,theta_nz))
    }

    plot(model_size,type = "l",ylab = "Model size",xlab = "iterations")


  }else if (type=="dist" & length(coef_val)==1){
    hist(ff$coef[,coef_val],main = "",xlab = paste0( expression(beta),"_",coef_val) )

  }else if (type=="dist" & length(coef_val)==2){
    hist(ff$coef_theta[,coef_val[1],coef_val[2]],main = "",xlab = paste0( expression(theta),"_",coef_val[1],"_",coef_val[2]) )

  }else if (type=="val" & length(coef_val)==1){
    plot(ff$pos_mpm_beta,main = "",xlab =  expression(beta) ,ylab = "coefficient")

  }else if (type=="val" & length(coef_val)==2){
    plot(c(ff$pos_mpm_theta),main = "",xlab = expression(theta),ylab = "coefficient" )

  }else if (type=="cont" ){
    xb<- ff$coef[,coef_val[1]]
    yt <- ff$coef_theta[,coef_val[1],coef_val[2]]
    s <- subplot(
      plot_ly(x = xb, type = "histogram",color = I("red")),
      plotly_empty(),
      plot_ly(x = xb, y = yt, type = "histogram2dcontour"),
      plot_ly(y = yt, type = "histogram",color = I("blue")),
      nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
      shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
    )
    fig <- layout(s, showlegend = FALSE)

    fig
    #hist(ff$coef[,coef_val],main = "",xlab = paste0( expression(beta),"_",coef_val) )

  }

}


