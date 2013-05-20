#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double sumC(NumericVector x) {
  double ans;
  for (int i = 0; i < x.size(); i++) {
    ans += x[i];
  }
  return ans;
}

// [[Rcpp::export]]
double sumC2(NumericVector x) {
  double ans = sum(x);
  return ans;
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
NumericVector pdistC2(double x, NumericVector ys) {
  NumericVector ans(ys.size());
  ans = sqrt(pow(x - ys, 2));  
  return ans;
}


// [[Rcpp::export]]
NumericVector rowSumC(NumericMatrix x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);  
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      ans[i] += x(i,j);
    }
  }  
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
LogicalVector biggerY(NumericVector x, double y) {
  LogicalVector ans(x.size());
  for (int i = 0; i < x.size(); i++) {
    ans[i] = x[i] > y;
  }
  return ans;
}

// [[Rcpp::export]]
LogicalVector biggerY2(NumericVector x, double y) {
  return x > y;
}


// [[Rcpp::export]]
NumericVector rollP3(NumericVector A, int start, int end) {
  NumericVector B = A[seq(start-1, end-1)] ;
  return B.sort() ;
}

// [[Rcpp::export]]
List scalar_missings() {
  int int_s = NA_INTEGER;
  String chr_s = NA_STRING;
  bool lgl_s = NA_LOGICAL;
  double num_s = NA_REAL;

  return List::create(int_s, chr_s, lgl_s, num_s);
}


// [[Rcpp::export]]
NumericVector biggerYsubset(NumericVector x, double y) {
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


// foo <- function(x, y) ifelse(x < y, x*x, -(y*y))

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


// [[Rcpp::export]]
NumericVector foov(NumericVector xx, NumericVector yy) {
  NumericVector res1 = ifelse(xx < yy, xx * xx, -(yy * yy));
}


double square( double x ){
  return x*x ;
}

// [[Rcpp::export]]
NumericVector applyC(NumericVector xx){
  return sapply(xx, square);
}

// fr <- function(n) {
//   if (n < 2) return(n)
//   return(fr(n-1) + fr(n-2))
// }

// [[Rcpp::export]]
int fC(int n) {
  if (n < 2) return(n);
  else return(fC(n-1) + fC(n-2));
}

// benchmark(fr(20), fC(20))[, 1:4]



// from http://gallery.rcpp.org/articles/vector-cumulative-sum/

#include <Rcpp.h>
#include <numeric> // for std::partial_sum

using namespace Rcpp;

// [[Rcpp::export]]
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
NumericVector cumsum2(NumericVector x) {
  NumericVector res(x.size());
  std::partial_sum(x.begin(), x.end(),
		   res.begin());
  return res;
}

// [[Rcpp::export]]
NumericVector cumsum3(NumericVector x){
  return cumsum(x); // compute + return result vector
}


// [[Rcpp::export]]
NumericVector callFunction(NumericVector x,
			   Function f) {
  NumericVector res = f(x);
  return res;
}

// callFunction(x, summary)

