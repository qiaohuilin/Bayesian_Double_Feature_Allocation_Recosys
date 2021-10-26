#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix readintb(NumericMatrix tb, NumericMatrix R, int dim1, int dim2, int sta, int rl){
    for(int i=sta; i<sta+rl; i++){
      int m= tb(i,0);
      int n= tb(i,1);
      R(m-1,n-1)=tb(i,2);
    }
  return(R);
}


NumericMatrix readintb1(NumericMatrix tb, NumericMatrix R, int dim1, int dim2, int sta, int rl){
  for(int i=sta; i<sta+rl; i++){
    int m= tb(i,0);
    int n= tb(i,1);
    R(m,n)=tb(i,2);
  }
  return(R);
}