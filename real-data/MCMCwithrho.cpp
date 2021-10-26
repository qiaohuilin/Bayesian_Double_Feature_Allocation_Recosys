#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
IntegerVector k_ij(NumericVector Ai, NumericVector Bj){
  IntegerVector idx= seq(0,Ai.size()-1);
  IntegerVector k=idx[Ai!=0 & Bj!=0];
  return(k);
}

// [[Rcpp::export]]
double mu_ij(NumericVector Ai, NumericVector Bj, NumericVector theta){
  IntegerVector k=k_ij(Ai,Bj);
  NumericVector thetak=theta[k];
  double mu=sum(thetak);
  return mu;
}

// [[Rcpp::export]]
double Likelihood(double y_ij,NumericVector Ai, NumericVector Bj, NumericVector theta, double tau, double rhoj){
  IntegerVector k=k_ij(Ai,Bj);
  NumericVector thetak=theta[k];
  double ll  = 0.0 ;
  if(y_ij==1){
    ll= R::pnorm(1,sum(thetak)+rhoj,tau, 1, 0);
  }
  if(y_ij>=2 & y_ij<=4){
    ll=R::pnorm(y_ij,sum(thetak)+rhoj,tau,1, 0) - R::pnorm((y_ij-1),sum(thetak)+rhoj,tau,1,0);
  }
  if(y_ij==5){
    ll= 1- R::pnorm(4,sum(thetak)+rhoj,tau, 1, 0);
  }
  return (ll);
}

// [[Rcpp::export]]
double Likelihood_oneuser(int i, int n, NumericMatrix R_obs, NumericMatrix A, NumericMatrix B, NumericVector theta, double tau, NumericVector rho){
  int id= i-1;
  double logL=0.0;
  double ll=0.0;
  for(int j=0; j<n; j++){
    if(R_obs(id,j)!=0){
      ll=Likelihood(R_obs(id,j),A(id,_),B(j,_),theta,tau,rho[j]);
      logL=logL + log(ll);
    }
  }
  return(exp(logL));
}

// [[Rcpp::export]]
double Likelihood_onemovie(int j, int m, NumericMatrix R_obs, NumericMatrix A, NumericMatrix B, NumericVector theta, double tau, NumericVector rho){
  int jd= j-1;
  double logL=0.0;
  double ll=0.0;
  for(int i=0; i<m; i++){
    if(R_obs(i,jd)!=0){
      ll=Likelihood(R_obs(i,jd),A(i,_),B(jd,_),theta,tau,rho[jd]);
      logL=logL + log(ll);
    }
  }
  return(exp(logL));
}

// [[Rcpp::export]]
double log_likelihood_all(int m,int n,NumericMatrix R_obs, NumericMatrix A, NumericMatrix B, NumericVector theta, double tau, NumericVector rho){
  double L= 0.0;
  double ll=0.0;
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      if(R_obs(i,j)!=0){
        ll=Likelihood(R_obs(i,j),A(i,_),B(j,_),theta,tau,rho[j]);
        L=L+log(ll);
      }
    }
  }
  return(L);
}

//[[Rcpp::export]]
double rbern1(double p){
  NumericVector y=Rcpp::rbinom(1,1,p);
  return(y[0]);
}

// [[Rcpp::export]]
NumericMatrix update_B(int K_tot, int m, int n, NumericMatrix R_obs, NumericMatrix A, NumericVector theta, double tau,NumericVector rho){
  NumericMatrix B(n, (K_tot+1));
  for(int j=0; j<n;j++){
    B(j,0)=1;
  }
  double p0=0.0;
  double p1=0.0;
  for(int j=0; j<n; j++){
    for(int k=1; k<(K_tot+1); k++){
      NumericMatrix B0temp=B;
      B0temp(j,k)=0;
      p0= 0.8 * Likelihood_onemovie(j+1,m,R_obs,A,B0temp,theta,tau,rho);
      
      NumericMatrix B1temp=B;
      B1temp(j,k)=1;
      p1= 0.2 * Likelihood_onemovie(j+1,m,R_obs,A,B1temp,theta,tau,rho);
      
      if((p0==0.0 & p1==0.0) | (R_IsNA(p0)&R_IsNA(p1))){
        p1=0.8;
        p0=0.2;
      }else{
        p1=p1/(p0+p1);
      }
      
      B(j,k)=rbern1(p1);
    }
  }
  
  return(B);
}

// [[Rcpp::export]]
List update_Z_tau(int m,int n, NumericMatrix R_obs, NumericMatrix A, NumericMatrix B, NumericVector theta, double tau, NumericVector rho){
  List ret;
  NumericMatrix Z(m,n);
  double err=0.0;
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      if(R_obs(i,j)!=0){
        IntegerVector k=k_ij(A(i,_),B(j,_));
        NumericVector thetak=theta[k];
        double mu=sum(thetak);
        if(R_obs(i,j)==1){
          Z(i,j)=r_truncnorm(mu+rho[j],tau,-pow(10,5),1);
        }
        if(R_obs(i,j)>=2 & R_obs(i,j)<=4){
          Z(i,j)=r_truncnorm(mu+rho[j],tau,R_obs(i,j)-1, R_obs(i,j));
        }
        if(R_obs(i,j)==5){
          Z(i,j)=r_truncnorm(mu+rho[j],tau,4, pow(10,5));
        }
        
        err=err+pow((Z(i,j)-mu-rho[j]),2);
      }
    }
  }
  double tau_update=sqrt(1/R::rgamma(5+1/2*sum(R_obs!=0), 1/(1+1/2*err)));
  ret["Z"]=Z;
  ret["err"]=err;
  ret["tau"]=tau_update;
  return ret;
}


// [[Rcpp::export]]
NumericVector update_theta(int K_tot, int m, int n, NumericMatrix R_obs, NumericMatrix A, NumericMatrix B, NumericMatrix Z, NumericVector theta0, double tau, double sig0, NumericVector rho){
  NumericVector theta=theta0;
  int nk=0;
  double mu=0;
  IntegerVector kij; 
  int k=1;
  while(k<K_tot){
    nk=0;
    mu=0;
    for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
        if(R_obs(i,j)!=0){
          kij=k_ij(A(i,_),B(j,_));
          if(std::find(kij.begin(), kij.end(), k) != kij.end() ){
            nk=nk+1;
            NumericVector thetakij=theta[kij];
            mu=mu+ (Z(i,j)-(sum(thetakij)-theta[k])-rho[j]);
          }
        }
      }
    }
    double t1=pow(sig0,2);
    double t2=pow(tau,2);
    theta[k]=R::rnorm((mu/t2)/(1/t1+nk/t2),sqrt(1/(1/t1+nk/t2)));
    k++;
  }
  return(theta);
}

// [[Rcpp::export]]
NumericVector update_rho(int m, int n, NumericMatrix R_obs, NumericMatrix A, NumericMatrix B, NumericMatrix Z, NumericVector theta, double tau,NumericVector rho0){
  NumericVector rho=rho0;
  IntegerVector kij;
  NumericVector thetak;
  double mu=0;
  int nj=0;
  for(int j=0;j<n;j++){
    mu=0;
    nj=0;
    for(int i=0; i<m;i++){
      if(R_obs(i,j)!=0){
        kij=k_ij(A(i,_),B(j,_));
        thetak=theta[kij];
        mu=mu+ (Z(i,j)-sum(thetak));
        nj=nj+1;
      }
    }
    double t1=pow(1,2);
    double t2=pow(tau,2);
    rho[j]= R::rnorm((mu/t2)/(1/t1+nj/t2),sqrt(1/(1/t1+nj/t2)));
  }
  
  return(rho);
}


// [[Rcpp::export]]
double tntest(double mu,double sigma,double a,double b){
  double r= r_truncnorm(mu,sigma,a,b);
  return(r);
}

double rbern(double p){
  NumericVector y=Rcpp::rbinom(1,1,p);
  return(y[0]);
}

double rnormtest(double sig0,double tau,double mu,int nk){
  double t1=pow(sig0,2);
  double t2=pow(tau,2);
  double th=R::rnorm((mu/t2)/(1/t1+nk/t2),sqrt(1/(1/t1+nk/t2)));
  return(th);
}

int isin(NumericVector kij){
  int y=0;
  if(std::find(kij.begin(), kij.end(), 0) != kij.end() ){
    y=1;
  }else{
    y=0;
  }
  return(y);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


