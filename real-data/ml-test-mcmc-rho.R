rm(list=ls())
library(MCMCpack)
library(invgamma)
library(truncnorm)
library(extraDistr)
library(EnvStats)

#Rcpp::sourceCpp("~/Desktop/Peter/simulations/MCMCforsim5_ord_cpp.cpp")
Rcpp::sourceCpp("Movielens/MCMCwithrho.cpp")

set.seed(123)
load("movielensrating.RData")

###for shard 1 ########
R_true=R1[1:400,]
m=150
n=200

R_obs=R_true
for(i in 1:m){
    R_obs[i,sample(which(R_true[i,]!=0),size=1)]=0
}
train_bi=1*(R_obs!=0)
test_bi=(1*(R_true!=0))*(1*(R_obs==0))
R_test=R_true*test_bi

#Initialize A, B, K using matrix factorization with a grid number of features
#in this section, lambda is regularization, alpha is step size.
unroll_Vecs <- function (params, Y, R, num_users, num_movies, num_features) {
  # Unrolls vector into X and Theta
  # Also calculates difference between preduction and actual 
  
  endIdx <- num_movies * num_features
  
  X     <- matrix(params[1:endIdx], nrow = num_movies, ncol = num_features)
  Theta <- matrix(params[(endIdx + 1): (endIdx + (num_users * num_features))], 
                  nrow = num_users, ncol = num_features)
  
  Y_dash     <-   (((X %*% t(Theta)) - Y) * R) # Prediction error
  
  return(list(X = X, Theta = Theta, Y_dash = Y_dash))
}

J_cost <-  function(params, Y, R, num_users, num_movies, num_features, lambda, alpha) {
  # Calculates the cost
  
  unrolled <- unroll_Vecs(params, Y, R, num_users, num_movies, num_features)
  X <- unrolled$X
  Theta <- unrolled$Theta
  Y_dash <- unrolled$Y_dash
  
  J <-  .5 * sum(   Y_dash ^2)  + lambda/2 * sum(Theta^2) + lambda/2 * sum(X^2)
  
  return (J)
}

grr <- function(params, Y, R, num_users, num_movies, num_features, lambda, alpha) {
  # Calculates the gradient step
  # Here lambda is the regularization parameter
  # Alpha is the step size

  unrolled <- unroll_Vecs(params, Y, R, num_users, num_movies, num_features)
  X <- unrolled$X
  Theta <- unrolled$Theta
  Y_dash <- unrolled$Y_dash
  
  X_grad     <- ((   Y_dash  %*% Theta) + lambda * X     )
  Theta_grad <- (( t(Y_dash) %*% X)     + lambda * Theta )
  
  grad = c(X_grad, Theta_grad)
  return(grad)
}

num_features_list=c(10,15,20,25,30,35)
cost <- NULL
for(i in 1:4){
    num_features=num_features_list[i]
    res <- optim(par = c(runif(m * num_features), runif(n * num_features)), # Random starting parameters
               fn = J_cost, gr = grr, 
               Y=t(R_obs), R=t(train_bi), 
               num_users=m, num_movies=n,num_features=num_features, 
               lambda=1, alpha = 0.2, 
               method = "L-BFGS-B", control=list(maxit=3000, trace=1))
    cost[i]=res$value
}
cost[2:4]-cost[1:3]

num_features=15
res <- optim(par = c(runif(m * num_features), runif(n * num_features)), # Random starting parameters
               fn = J_cost, gr = grr, 
               Y=t(R_obs), R=t(train_bi), 
               num_users=m, num_movies=n,num_features=num_features, 
               lambda=1, alpha = 0.2, 
               method = "L-BFGS-B", control=list(maxit=3000, trace=1))
U=matrix(res$par[1:m * num_features],m,num_features)
V=matrix(res$par[1:n * num_features],n,num_features)

A=1*(U>0.7)
B=1*(V>0.7)
K_tot=15
A=cbind(rep(1,m),A)
B=cbind(rep(1,n),B)

#Initialize
lambda=3 #mass parameter for IBP, noted as gamma in T.B review paper, or alpha in the original IBP paper
sig0=2
tau=sqrt(rinvgamma(1,5,1))


theta=c(2.5,rnorm(n=K_tot,mean=0,sd=sig0))


#add rho from prior N(0,1), initialize at rho_j=mean(R[,j])-2.5
rho_init=sapply(1:n, function(x) mean(R_obs[,x][R_obs[,x]!=0]))-2.5
rho=rho_init

#train theta, tau rho first for a better initialization, at current A,B,K
init_iter=0
while(init_iter<1000){
  ret=update_Z_tau(m,n,R_obs,A,B,theta,tau,rho)
  Z=ret$Z
  tau=sqrt(rinvgamma(1,5+1/2*sum(R_obs!=0),1+1/2*ret$err))
  theta=update_theta(K_tot, m, n, R_obs, A, B,Z, theta, tau, sig0,rho)
  rho=update_rho(m, n, R_obs, A, B,Z, theta, tau, rho)
  init_iter=init_iter+1
}


#MCMC update
# functions of mu_ij, Likelihood, Likelihood_oneuser, Likelihood_onemovies,
# log_likelihood_all are from cpp

tot_iter=5000
iter=1
tt=1
theta_store=vector("list",(tot_iter-805)/5+1)
tau_store=vector("list",(tot_iter-805)/5+1)
A_store=vector("list", (tot_iter-805)/5+1)
B_store=vector("list", (tot_iter-805)/5+1)
k_store <- NULL
rho_store=matrix(0,nrow=(tot_iter-805)/5+1,ncol=n)
log_likelihood_train_store <- NULL
log_likelihood_test_store <- NULL
R_pred=array(0,dim=c(m,n,5))
mu_pred_mat=matrix(0,m,n)



while(iter<tot_iter+1){
  #at iteration iter
  #update A(t+1)
  for (i in 1:m){
    for (k in 2: dim(A)[2]){ # +1 if for the structure of while, so dim(A)[2]th column still gets run
      if(sum(A[-i,k])!=0){
        A0temp=A
        A0temp[i,k]=0
        p0= (1-sum(A[-i,k])/m) * Likelihood_oneuser(i,n,R_obs,A0temp,B,theta,tau,rho)
        
        A1temp=A
        A1temp[i,k]=1
        p1= sum(A[-i,k])/m * Likelihood_oneuser(i,n,R_obs,A1temp,B,theta,tau,rho)
        
        if((p0==0 & p1==0) | (is.na(p0) & is.na(p1))){
          p1=1-sum(A[-i,k])/m
          p0=sum(A[-i,k])/m
        }else{
          p1=p1/(p0+p1)
        }
        A[i,k]= rbern(1,prob=p1)
      }
    }
    
    Atemp=A
    Btemp=B
    thetatemp=theta
    K_tottemp=K_tot
    for (k in  2: dim(A)[2]){ # +1 if for the structure of while, so dim(A)[2]th column still gets run
      de=0
      if(all(A[-i,k]==0)){
        Atemp=as.matrix(A[,-k])
        Btemp=as.matrix(B[,-k])
        thetatemp=theta[-k]
        K_tottemp=K_tot-1
        de=de+1
      }
    }
    
    k_star=rpois(n=1,lambda/m)
    K_tottemp=K_tottemp+k_star
    
    if(k_star>0){
      Btemp=cbind(Btemp,matrix(rbern(n=n*k_star,prob=0.5),nrow=n,ncol=k_star))
      Atemp=cbind(Atemp,matrix(0,nrow=m,ncol=k_star))
      Atemp[i,(K_tottemp+1-k_star+1):(K_tottemp+1)]=1   #original dim is K_tot +1 the intercept 
      thetatemp=c(thetatemp,rnorm(n=k_star,mean=0,sd=2))
    }
    acc_rate=Likelihood_oneuser(i,n,R_obs,Atemp,Btemp,thetatemp,tau,rho)/Likelihood_oneuser(i,n,R_obs,A,B,theta,tau,rho)
    if(is.na(acc_rate)){
      acc_rate=1
    }
    acc_rate=min(1,acc_rate)
    
    u=runif(1,0,1)
    if(u < acc_rate){
      A=Atemp
      B=Btemp
      theta=thetatemp
      K_tot=K_tottemp
    }
    
    #print(c(de,k_star))
  }
  
  #update B
  # use Rcpp function
  B=update_B(K_tot, m, n,  R_obs, A, theta, tau,rho)
  
  # update Z (truncated normal given y) and tau^2 invgamma
  # use Rcpp function
  ret=update_Z_tau(m,n,R_obs,A,B,theta,tau,rho)
  Z=ret$Z
  tau=sqrt(rinvgamma(1,5+1/2*sum(R_obs!=0),1+1/2*ret$err))
  
  # update theta
  # use Rcpp function
  theta=update_theta(K_tot, m, n, R_obs, A, B,Z, theta, tau, sig0,rho)
  
  #update rho
  rho=update_rho(m, n, R_obs, A, B,Z, theta, tau, rho)
  
  
  log_likelihood_train_store[iter]=log_likelihood_all(m,n,R_obs,A,B,theta,tau,rho)
  #log_likelihood_test_store[iter]=log_likelihood_all(m,n,R_mis,A,B,theta,tau)
  
  print(c(iter,K_tot))
  
  #prediction
  if(iter>800 & iter%%5==0){
    for (i in 1:m){
      for (j in 1:n){   
        mu_pred_mat[i,j]=mu_pred_mat[i,j]+mu_ij(A[i,],B[j,],theta) 
        R_pred[i,j,]= R_pred[i,j,]+ sapply(1:5, function(x) Likelihood(x, A[i,],B[j,],theta,tau,rho[j]))
      }
    }
    A_store[[tt]] <- A
    B_store[[tt]] <- B
    k_store[[tt]] <- K_tot
    theta_store[[tt]]  <- theta
    rho_store[tt,]<- rho
    tau_store[[tt]] <- tau
    print(rho_store[tt,])
    tt=tt+1
  }
  
  iter=iter+1
}


R_pred_vote= matrix(0,m,n)
for(i in 1:m){
  for(j in 1:n){
    R_pred_vote[i,j]=which.max(R_pred[i,j,])
  }
}


tr=R_pred_vote * train_bi
ts=R_pred_vote * test_bi
sum(tr[tr !=0]==R_obs[R_obs!=0])/sum(train_bi)
sum(ts[ts !=0]==R_test[R_test!=0])/sum(test_bi)


#compare pairs
pair_true=NULL
idx1=NULL
idx2=NULL
for(i in 1:m){
    idx1[i]=which(test_bi[i,]!=0)
    x1=R_true[i,idx1[i]]
    idx2[i]=sample(which(R_obs[i,]!=0),size=1)
    x2=R_true[i,idx2[i]]
    if(x1>x2){
        pair_true[i]=1
    }else{
        pair_true[i]=0
    }  
}

pair_pred=NULL
pair_pred2=NULL
for(i in 1:m){
    x1=R_pred_vote[i,idx1[i]]
    x2=R_pred_vote[i,idx2[i]]
    x22=R_true[i,idx2[i]]
    if(x1>x2){
        pair_pred[i]=1
    }else{
        pair_pred[i]=0
    }
    if(x1>x22){
        pair_pred2[i]=1
    }else{
        pair_pred2[i]=0
    }
}

1-sum((pair_true-pair_pred)^2)/length(pair_true)
1-sum((pair_true-pair_pred2)^2)/length(pair_true)


save.image("ml_withrho_shard1_m=150_n=200_withmu_rcpp_withAB.RData")

