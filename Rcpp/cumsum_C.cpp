#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;
//[[Rcpp::export]]

double cumsum_C(NumericVector x){
  return std::accumulate(x.begin(),x.end(),0.0);
}
