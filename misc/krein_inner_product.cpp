
#include <Rcpp.h>
using namespace Rcpp;




// [[Rcpp::export]]
NumericMatrix kreinProduct(NumericMatrix X, NumericMatrix Y,int p ) {
  //function to calculate pairwise krein inner product for
  //two matrices
  int n = X.nrow();
  int d = X.ncol();
  NumericMatrix toReturn(n,n);
  
  //iterate through all of X
  for (int i = 0; i < n; i ++) {
    //iterate through all of Y:
    for(int j=0;j<=i;j++) {
      //iterate through all of the dimensions
      for (int k = 0; k < d;k++) {
        if (k < p) {
          toReturn(i,j) += X(i,k)*Y(j,k);
        } else {
          toReturn(i,j) -= X(i,k)*Y(j,k);
        }
        
      }
      
      if(i != j) {
        toReturn(j,i) = toReturn(i,j); //symmetric
      }
      
      
    }
  }
  
  return toReturn;
  
}


