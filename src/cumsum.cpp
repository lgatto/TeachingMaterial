#include <Rcpp.h>
#include <numeric> // for std::partial_sum

using namespace Rcpp;

// from http://gallery.rcpp.org/articles/vector-cumulative-sum/
NumericVector cumsum1(NumericVector x){
  // initialize an accumulator variable
  double acc = 0;
  // initialize the result vector
  NumericVector res(x.size());
  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  return res;
}

// [[Rcpp::export]]
NumericVector cumsum2(NumericVector x){
  return cumsum(x); // compute + return result
}


// [[Rcpp::export]]
NumericVector cumsum3(NumericVector x) {
  NumericVector res(x.size());
  std::partial_sum(x.begin(), x.end(), res.begin());
  return res;
}
