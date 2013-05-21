#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector modifyX(IntegerVector x) {
  int n = x.size();
  for (int i = 0; i < n; i++) {
    x[i] = i;
  }
  return x;
}

// [[Rcpp::export]]
IntegerVector nomodifyX(IntegerVector x) {
  IntegerVector y = clone(x); 
  int n = x.size();
  for (int i = 0; i < n; i++) {
    y[i] = i;
  }
  return y;
}

