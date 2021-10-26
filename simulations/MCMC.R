rm(list=ls())
library(MCMCpack)
library(invgamma)
library(truncnorm)
library(extraDistr)
library(EnvStats)

Rcpp::sourceCpp("MCMC_cpp.cpp")

set.seed(123)
load("trueR5_lam=3_b=0.2_m=100_n=150.RData")
#load("trueR5_lam=3_b=0.5_m=50.RData")

R_obs=R_obs
R_mis=R*(1-obs_bi)

m=100
n=150

#Initialize at the true_value
lambda=3 #mass parameter for IBP, noted as gamma in T.B review paper, or alpha in the original IBP paper
#theta =1 concentration parameter in T.B review paper, ignored as 1 in other papers
A=A_true
B=B_true
K_tot=K_tot # does not include intercept! So dim(A)=K_tot+1

sig0=2
theta=c(2.5,rnorm(n=K_tot,mean=0,sd=sig0))
tau=sqrt(rinvgamma(1,5,1))

#MCMC update
# functions of mu_ij, Likelihood, Likelihood_oneuser, Likelihood_onemovies,
# log_likelihood_all are from cpp

tot_iter=8000
iter=1
tt=1
theta_store=vector("list",(tot_iter-805)/5+1)
tau_store=vector("list",(tot_iter-805)/5+1)
A_store=vector("list", (tot_iter-805)/5+1)
B_store=vector("list", (tot_iter-805)/5+1)
k_store <- NULL
log_likelihood_train_store <- NULL
log_likelihood_test_store <- NULL
R_pred=array(0,dim=c(m,n,5))
mu_pred_mat=matrix(0,m,n)

mylist <- vector("list",tot_iter)

while(iter<tot_iter+1){
  #at iteration iter
  #update A(t+1)
  for (i in 1:m){
    for (k in 2: dim(A)[2]){ # +1 if for the structure of while, so dim(A)[2]th column still gets run
      if(sum(A[-i,k])!=0){
        A0temp=A
        A0temp[i,k]=0
        p0= (1-sum(A[-i,k])/m) * Likelihood_oneuser(i,n,R_obs,A0temp,B,theta,tau)
        
        A1temp=A
        A1temp[i,k]=1
        p1= sum(A[-i,k])/m * Likelihood_oneuser(i,n,R_obs,A1temp,B,theta,tau)
        
        if(p0==0 & p1==0){
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
    acc_rate=Likelihood_oneuser(i,n,R_obs,Atemp,Btemp,thetatemp,tau)/Likelihood_oneuser(i,n,R_obs,A,B,theta,tau)
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
  B=update_B(K_tot, m, n,  R_obs, A, theta, tau)
  
  # update Z (truncated normal given y) and tau^2 invgamma
  # use Rcpp function
  ret=update_Z_tau(m,n,R_obs,A,B,theta,tau)
  Z=ret$Z
  tau=sqrt(rinvgamma(1,5+1/2*sum(R_obs!=0),1+1/2*ret$err))
  
  # update theta
  # use Rcpp function
  theta=update_theta(K_tot, m, n, R_obs, A, B,Z, theta, tau, sig0)
  
  
  theta_store[[iter]]  <- theta
  tau_store[[iter]] <- tau
  log_likelihood_train_store[iter]=log_likelihood_all(m,n,R_obs,A,B,theta,tau)
  log_likelihood_test_store[iter]=log_likelihood_all(m,n,R_mis,A,B,theta,tau)
  
  print(c(iter,K_tot))
  
  #prediction
  if(iter>800 & iter%%5==0){
    for (i in 1:m){
      for (j in 1:n){   
        mu_pred_mat[i,j]=mu_pred_mat[i,j]+mu_ij(A[i,],B[j,],theta) 
        R_pred[i,j,]= R_pred[i,j,]+ sapply(1:5, function(x) Likelihood(x, A[i,],B[j,],theta,tau))
      }
    }
    A_store[[tt]] <- A
    B_store[[tt]] <- B
    k_store[[tt]] <- K_tot
    theta_store[[tt]]  <- theta
    tau_store[[tt]] <- tau
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

tr=R_pred_vote * obs_bi
ts=R_pred_vote * (1-obs_bi)
sum(tr[tr !=0]==R_obs[R_obs!=0])/sum(obs_bi)
sum(ts[ts !=0]==R_mis[R_mis!=0])/sum(1-obs_bi)

save.image("sim5_lam=3_b=0.5_m=100_n=150_withmu_rcpp.RData")

#classification
sum(R_pred_vote==R)/(m*n)
#regression mse  (the book did use this residual)
sum((R_pred_vote-R)^2)/(m*n)
sum((R_pred_vote*obs_bi-R*obs_bi)^2)/(sum(obs_bi))
sum((R_pred_vote*(1-obs_bi)-R*(1-obs_bi))^2)/(sum(1-obs_bi))
#wilcoxon signed rank test
wilcox.test(R,R_pred_vote,mu=0,alt="two.sided",paired=T, conf.int=T,conf.level=0.9)
#wilcoxon rank sum test
wilcox.test(R,R_pred_vote,mu=0,alt="two.sided",paired=F, conf.int=T,conf.level=0.9)
#friedman test and kendall's W
friedman.test(R,R_pred_vote)
library(irr)
kendall(cbind(as.vector(R),as.vector(R_pred_vote)), correct = T)
kruskal.test(list(as.vector(R),as.vector(R_pred_vote)))
#residual deviance?? I don't seem to have a constant df q

#residual of zs and mu_ij
#mean((Z_true-mu_pred_mat/((iter-1-800)/5+1))^2)

#prob vec and euclidean dist
R_true_prob=array(0,dim=c(m,n,5))
R_pred_prob=R_pred/((iter-1-800)/5+1)
mu_pred_mat=mu_pred_mat/((iter-1-800)/5+1)
p0 = matrix(0,m,n)
pr = matrix(0,m,n)
D = matrix(0,m,n)

for(i in 1:m){
  for(j in 1:n){
    R_true_prob[i,j,]=sapply(1:5,function(x) Likelihood(x,A_true[i,],B_true[j,],theta_true,tau_true))
    p0[i,j]=R_true_prob[i,j,R[i,j]]
    pr[i,j]=R_pred_prob[i,j,R[i,j]]
    D[i,j]=sqrt(sum((R_true_prob[i,j,]-R_pred_prob[i,j,])^2)/5)
  }
}


mean(p0)
mean(pr)
mean(D)
sum(p0*obs_bi)/sum(obs_bi)
sum(pr*obs_bi)/sum(obs_bi)
sum(D*obs_bi)/sum(obs_bi)

p0_obs=p0*obs_bi
p0_obs=p0_obs[p0_obs!=0]

pr_obs=pr*obs_bi
pr_obs=pr_obs[pr_obs!=0]

D_obs=D*obs_bi
D_obs=D_obs[D!=0]

hist(p0_obs,prob=TRUE)
hist(pr_obs,prob=T)
hist(D_obs,prob=T)

R_true_vote= matrix(0,m,n)
c=0
for(i in 1:m){
  for(j in 1:n){
    if(R_obs[i,j]!=0){
      if(R_obs[i,j]==which.max(R_true_prob[i,j,])){c=c+1}
    }
  }
}
c/sum(obs_bi)


c=0
for(i in 1:m){
  for(j in 1:n){
    if(R_obs[i,j]!=0){
      if(R_obs[i,j]==which.max(R_pred_prob[i,j,])){c=c+1}
    }
  }
}
c/sum(obs_bi)

#if I use jensen=shannon divergence
#here I use log2 base so that they 0<jsd<1
H <- function(x){
  if(any(x==0)){
    return(-sum(x[x!=0]*log2(x[x!=0])))
  }else{
    return(-sum(x*log2(x)))
  }
}
jsd <- function(x,y){
  return(H(0.5*x+0.5*y)-0.5*H(x)-0.5*H(y))
}
jsd(R_true_prob[i,j,],R_pred_prob[i,j,])

#compare to the library
library(philentropy)
JSD(rbind(R_true_prob[i,j,],R_pred_prob[i,j,]),unit="log2")
JSD(rbind(R_true_prob[i,j,],R_pred_prob[i,j,]),unit="log")

#so for all i,j
jsd_mat=matrix(0,m,n)
for(i in 1:m){
  for(j in 1:n){
    jsd_mat[i,j]=jsd(R_true_prob[i,j,],R_pred_prob[i,j,])
  }
}
hist(jsd_mat)
hist(jsd_mat[jsd_mat<0.1],breaks=10,prob=T)

#rank and pairwise preference
rank_pred=matrix(0,m,n)
rank_true=matrix(0,m,n)
for(i in 1:m){
  rank_pred[i,]=n+1-rank(mu_pred_mat[i,])
  rank_true[i,]=n+1-rank(mu_true[i,])
}

pref_testpair_accuracy_predrate=NULL
pref_testpair_accuracy_predrank=NULL
for(i in 1:m){
  for(j in 1:n){
    if(obs_bi[i,j]==0){
      pref_pred_ij=sapply(rank_pred[i,-j], function(x) if(x<rank_pred[i,j]){return(0)}else{return(1)})
      pref_pred_ij2=sapply(R_pred_vote[i,-j], function(x) if(x>R_pred_vote[i,j]){return(0)}else{return(1)})
      pref_true_ij=sapply(R_true[i,-j], function(x) if(x>R_true[i,j]){return(0)}else{return(1)})
      pref_testpair_accuracy_predrank=c(pref_testpair_accuracy_predrank,sum(pref_pred_ij == pref_true_ij)/length(pref_pred_ij))
      pref_testpair_accuracy_predrate=c(pref_testpair_accuracy_predrate,sum(pref_pred_ij2 == pref_true_ij)/length(pref_pred_ij))
    }
  }
}

ftd=NULL
for(i in 1:m){
  ftd[i]=sum(abs(rank_pred[i,]-rank_true[i,]))
}
mean(ftd)


# if run mallows model on it (vitelli 2018)
library(BayesMallows)
library(dplyr)

#this is way too slow
pair_comp=tribble(~assessor, ~bottom_item, ~top_item,
                  0, 0, 0)
for(i in 1:m){
  print(i)
  for(j in which(obs_bi[i,]!=0)){
    for(k in which(obs_bi[i,]!=0)){
      if(j !=k){
        if(R_obs[i,j]>R_obs[i,k]){
          pair_comp=rbind(pair_comp,c(i,k,j))
        }else{
          pair_comp=rbind(pair_comp,c(i,j,k))
        }
      }
    }
  }
}

#let's try move to Rcpp
Rcpp::sourceCpp("~/Desktop/Peter/simulations/test_cpp.cpp")

pair_comp=pair(R_obs)
pair_comp=pair_comp[rowSums(pair_comp)!=0,]

pair_comp=data.frame(pair_comp)
colnames(pair_comp) <- c("assessor", "bottom_item", "top_item")

pair_tc=generate_transitive_closure(pair_comp, cl = NULL)

epf=estimate_partition_function(method = "importance_sampling", alpha_vector=seq(0,1,0.5),n_items=150, nmc=1000,metric="footrule",degree=2)
mallows_res=compute_mallows(rankings = NULL, preferences = pair_comp, metric = "footrule", error_model="bernoulli",nmc=50000,alpha_init = 1,alpha_jump=5,logz_estimate=epf,save_aug = TRUE)

assess_convergence(mallows_res)
#plot rho
plot(mallows_res,burnin = 1000,parameter = "rho")

#marginal Rtilde
assess_convergence(mallows_res, parameter = "Rtilde",items = c(2, 4), assessors = c(1, 2))

#mean posterior rank for assessor 1 item 1
burnin=1000
item_ind=levels(mallows_res$augmented_data$item)
rankij=mallows_res$augmented_data %>% filter(assessor==1) %>% filter(item %in% item_ind[1])
rankij=rankij$value[(burnin+1):length(rankij)]
rankij_postmean=mean(rankij)

#get mean rank for everyone
sm=mallows_res$augmented_data %>% group_by(assessor, item) %>%
  summarize(mean_rank = mean(value, na.rm = TRUE))
rankij_postmean_mat=matrix(sm$mean_rank,nrow=m,ncol=n,byrow = T)

#pairwise compare
pref_testpair_accuracy_predrank_mallows=NULL
for(i in 1:m){
  for(j in 1:n){
    if(obs_bi[i,j]==0){
      pref_pred_ij_ml=sapply(rankij_postmean_mat[i,-j], function(x) if(x<rankij_postmean_mat[i,j]){return(0)}else{return(1)})
      pref_true_ij=sapply(R_true[i,-j], function(x) if(x>R_true[i,j]){return(0)}else{return(1)})
      pref_testpair_accuracy_predrank_mallows=c(pref_testpair_accuracy_predrank_mallows,sum(pref_pred_ij_ml == pref_true_ij)/length(pref_pred_ij_ml))
    }
  }
}
mean(pref_testpair_accuracy_predrank_mallows)



