rm(list=ls())
library(MCMCpack)
library(invgamma)
library(extraDistr)

set.seed(101)
#simulate a 50 people 6 movie community. 
m=100
n=150

#simulate A mat
lambda=3 #mass parameter for IBP, noted as gamma in T.B review paper, or alpha in the original IBP paper
#theta =1 concentration parameter in T.B review paper, ignored as 1 in other papers
#propose A=50, never gets there
A=matrix(0,nrow=m,ncol=50)
K <- NULL
K[1]=rpois(n=1,lambda)
K_tot=K[1]
A[1,1:K[1]]=1
for(i in 2:m){
    for (j in 1:K_tot){
        A[i,j]=rbern(1,prob=sum(A[1:(i-1),j])/i)
    }
    K[i]=rpois(n=1,lambda/i)
    if(K[i]>0){
        A[i,(K_tot+1):(K_tot+K[i])]=1
        K_tot=K_tot+K[i]
    }
}
A=A[,1:K_tot]

#simulate B mat
B=matrix(0,nrow=n,ncol=K_tot)
for(i in 1:n){
    for(j in 1:K_tot){
        B[i,j]=rbern(1,prob=0.2)
    }
}

#give one column of all 1 to A and B  (the intercept)
A=cbind(rep(1,m),A)
B=cbind(rep(1,n),B)

tau=0.25
theta=c(2.5,rnorm(n=K_tot,mean=0,sd=2))

mu_true=matrix(0,nrow=m,ncol=n)
Z=matrix(0,nrow=m,ncol=n)
R=matrix(0,nrow=m,ncol=n)
for(i in 1:m){
    for(j in 1:n){
        k=which(A[i,]!=0 & B[j,]!=0)
        mu_true[i,j]=sum(theta[k])
        Z[i,j]=rnorm(1, mean=sum(theta[k]), sd=tau)
        if(Z[i,j]<=1){
            R[i,j]=1
        }
        if(Z[i,j]>1 & Z[i,j]<=2){
            R[i,j]=2
        }
        if(Z[i,j]>2 & Z[i,j]<=3){
            R[i,j]=3
        }
        if(Z[i,j]>3 & Z[i,j]<=4){
            R[i,j]=4
        }
        if(Z[i,j]> 4){
            R[i,j]= 5
        }
    }
}
table(R)

R_true=R
obs_bi=matrix(rbern(m*n,0.7),nrow=m,ncol=n)
R_obs=R*obs_bi

R_mis=R*(1-obs_bi)
A_true=A
B_true=B
tau_true=tau
theta_true=theta
Z_true=Z

save.image("trueR5_lam=3_b=0.2_m=100_n=150_tau=025.RData")


