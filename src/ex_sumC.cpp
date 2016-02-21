#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sumC(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}

/*** R
# Comparing with R
(x <- c(1, 3, rnorm(10)))
sumC(x)
sum(x)
*/
