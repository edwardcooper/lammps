#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
double mean_C(NumericVector x){
  int n =x.size();
  double total=0;
  for (int i=0;i<n;i++){
    total+=x[i];
  }
  return total/n;
}

/*** R
x=rnorm(1e7)
library(microbenchmark)
microbenchmark( mean(x),mean_C(x), times=1000)
*/
