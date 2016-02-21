#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sumC2(NumericVector x) {
  double ans = sum(x);
  return ans;
}


// [[Rcpp::export]]
NumericVector rowSumC2(NumericMatrix x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);  
  for (int i = 0; i < nr; i++) {
    ans[i] = sum(x(i,_));
  }  
  return ans;
}


// [[Rcpp::export]]
NumericVector pdistC2(double x, NumericVector ys) {
  NumericVector ans(ys.size());
  ans = sqrt(pow(x - ys, 2));  
  return ans;
}

// [[Rcpp::export]]
LogicalVector lgl_biggerYC2(NumericVector x, double y) {
  return x > y;
}


// [[Rcpp::export]]
NumericVector foov(NumericVector xx, NumericVector yy) {
  NumericVector res1 = ifelse(xx < yy, xx * xx, -(yy * yy));
}

