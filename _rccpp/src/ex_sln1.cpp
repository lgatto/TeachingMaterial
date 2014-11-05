#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int fibC(int n) {
  if (n == 0) return 0;
  if (n == 1) return 1;
  return fibC(n - 1) + fibC(n - 2);
}

// [[Rcpp::export]]
NumericVector pdistC(double x, NumericVector ys) {
  NumericVector ans(ys.size());
  for (int i = 0; i < ys.size(); i++) {
    ans[i] = sqrt(pow(x - ys[i], 2));
  }
  return ans;
}

// [[Rcpp::export]]
NumericVector setyC(NumericVector x, double y) {
  NumericVector ans(x.size());
  for (int i = 0; i < x.size(); i++) {
    if (x[i] > 0) {
      ans[i] = y;
    } else if (x[i] < 0) {
      ans[i] = -y;    
    }
  }
  return ans;
}

// [[Rcpp::export]]
LogicalVector lgl_biggerYC(NumericVector x, double y) {
  LogicalVector ans(x.size());
  for (int i = 0; i < x.size(); i++) {
    ans[i] = x[i] > y;
  }
  return ans;
}

// [[Rcpp::export]]
NumericVector biggerYC(NumericVector x, double y) {
  NumericVector tmp(x.size()), ans;
  int k = 0;
  for (int i = 0; i < x.size(); i++) {
    if (x[i] > y) { 
      tmp[k] = x[i];
      k++;
    }
  }
  if (k > 0) {
    ans = tmp[seq(0, k-1)];
  } else {
    ans[0] = NA_REAL;
  }
  return ans;
}


// Example from Dirk's slides

// [[Rcpp::export]]
NumericVector foo(NumericVector xx, NumericVector yy) {
  int n = xx.size();
  NumericVector res1( n );
  double x_ = 0.0, y_ = 0.0;
  for (int i = 0; i < n; i++) {
    x_ = xx[i];
    y_ = yy[i];
    if (R_IsNA(x_) || R_IsNA(y_)) {
      res1[i] = NA_REAL;
    } else if (x_ < y_) {
      res1[i] = x_ * x_;
    } else {
      res1[i] = -(y_ * y_);
    }
  }
  return(res1);
}

