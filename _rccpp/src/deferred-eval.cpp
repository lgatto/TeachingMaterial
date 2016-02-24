#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cond_sum_cpp(NumericVector x, NumericVector y, 
                    LogicalVector z) {
  double sum = 0;
  int n = x.length();

  for(int i = 0; i < n; i++) {
    if (!z[i]) continue;
    sum += x[i] + y[i];
  }
  return sum;
}

// [[Rcpp::export]]
double cond_sum_cpp2(NumericVector x, NumericVector y, 
		     LogicalVector z) {
  double sum = 0;
  int n = x.length();

  for(int i = 0; i < n; i++) {
    if (z[i]) sum += x[i] + y[i];    
  }
  return sum;
}
