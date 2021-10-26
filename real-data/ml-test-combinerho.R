rho_all=NULL
load("ml_withrho_shard1_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rho_store
load("ml_withrho_shard2_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard3_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard4_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard5_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard6_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard7_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard8_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard9_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)
load("ml_withrho_shard10_m=150_n=200_withmu_rcpp_withAB.RData")
rho_all=rbind(rho_all,rho_store)

hist(rho_all[,1])
hist(rho_all[,2])
hist(rho_all[,3])
....

rho_avg=colMeans(rho_all)

load("ml_withrho_shard1_m=150_n=200_withmu_rcpp_withAB.RData")
Rcpp::sourceCpp("~/Desktop/Peter/Movielens/MCMCwithrho.cpp")

R_pred_new=array(0,dim=c(m,n,5))
for(tt in 1:dim(rho_store)[1]){
  A=A_store[[tt]]
  B=B_store[[tt]]
  theta=theta_store[[tt]]
  tau=tau_store[[tt]]
  for (i in 1:m){
    for (j in 1:n){  
      R_pred_new[i,j,]= R_pred_new[i,j,]+ sapply(1:5, function(x) Likelihood(x, A[i,],B[j,],theta,tau,rho_avg[j]))
    }}
  print(tt)
}

R_pred_vote_new= matrix(0,m,n)
for(i in 1:m){
  for(j in 1:n){
    R_pred_vote_new[i,j]=which.max(R_pred_new[i,j,])
  }
}


tr=R_pred_vote_new * train_bi
ts=R_pred_vote_new * test_bi
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
  x1=R_pred_vote_new[i,idx1[i]]
  x2=R_pred_vote_new[i,idx2[i]]
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

save(R_pred_vote_new,file="ml-shard3-glpred-5shards.RData")

mat=matrix(0,2,5)
for(file in 1:10){
  filename=sprintf("ml_withrho_shard%i_m=150_n=200_withmu_rcpp_withAB.RData",file)
  load(filename)
  mat[1,file]=1-sum((pair_true-pair_pred2)^2)/length(pair_true)
}

for(file in 1:10){
  filename1=sprintf("ml_withrho_shard%i_m=150_n=200_withmu_rcpp_withAB.RData",file)
  filename2=sprintf("ml-shard%i-glpred-5shards.RData",file)
  load(filename1)
  load(filename2)
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
    x1=R_pred_vote_new[i,idx1[i]]
    x2=R_pred_vote_new[i,idx2[i]]
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
  
  mat[2,file]=1-sum((pair_true-pair_pred2)^2)/length(pair_true)
}

mat
save(mat,file="ml-combinerhohist-5shard.RData")
########################################################
rm(list=ls())
n=200
mu_rho=matrix(0,15,n)
prec_rho=matrix(0,15,n)

mat=matrix(0,2,15)

for(file in 1:15){
  filename=sprintf("ml_withrho_shard%i_m=400_n=200_withmu_rcpp_withAB.RData",file)
  load(filename)
  mat[1,file]=1-sum((pair_true-pair_pred2)^2)/length(pair_true)
  mu_rho[file,]=colMeans(rho_store)
  prec_rho[file,]=sapply(1:200, function(x) 1/var(rho_store[,x]))
}


mu_rho_all=colMeans(mu_rho)
prec_rho_all=colSums(prec_rho)

load("ml_withrho_shard12_m=400_n=200_withmu_rcpp_withAB.RData")
Rcpp::sourceCpp("Movielens/MCMCwithrho.cpp")

R_pred_new=array(0,dim=c(m,n,5))
for(tt in 1:dim(rho_store)[1]){
  A=A_store[[tt]]
  B=B_store[[tt]]
  theta=theta_store[[tt]]
  tau=tau_store[[tt]]
  for (i in 1:m){
    for (j in 1:n){  
      rhoj=rnorm(1,mu_rho_all[j],sd = sqrt(1/(prec_rho_all[j])))
      R_pred_new[i,j,]= R_pred_new[i,j,]+ sapply(1:5, function(x) Likelihood(x, A[i,],B[j,],theta,tau,rhoj))
    }}
  print(tt)
}

R_pred_vote_new= matrix(0,m,n)
for(i in 1:m){
  for(j in 1:n){
    R_pred_vote_new[i,j]=which.max(R_pred_new[i,j,])
  }
}


tr=R_pred_vote_new * train_bi
ts=R_pred_vote_new * test_bi
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
  x1=R_pred_vote_new[i,idx1[i]]
  x2=R_pred_vote_new[i,idx2[i]]
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
ts=1-sum((pair_true-pair_pred2)^2)/length(pair_true)

pred_new=list(R_pred_new=R_pred_new,ts=ts)

save(pred_new,file="ml-shard10-glpred-samplerho.RData")

mat=matrix(0,2,10)

for(file in 1:10){
  filename=sprintf("ml_withrho_shard%i_m=150_n=200_withmu_rcpp_withAB.RData",file)
  load(filename)
  mat[1,file]=1-sum((pair_true-pair_pred2)^2)/length(pair_true)
}


for(file in 1:10){
  filename2=sprintf("ml-shard%i-glpred-10shards-samplerho.RData",file)
  load(filename2)
  mat[2,file]=pred_new$ts
}

save(mat,file="ml-combinerhohist-10shard.RData")

