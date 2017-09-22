#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
List lapply_C(List input,Function f){
  int n =input.size();
  List out(n);
  for (int i=0;i<n;i++){
  out[i]=f(input[i]);
  }
  return out;
}

/*** R
library(microbenchmark)
long_list=list(1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6)
microbenchmark(lapply(long_list,mean),lapply_C(long_list,mean),times=1000)
*/
