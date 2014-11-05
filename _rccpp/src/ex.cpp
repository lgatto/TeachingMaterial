#include <Rcpp.h>

using namespace Rcpp;

////////////////////////
//     Exercise 1     //
////////////////////////

// [[Rcpp::export]]
double sumC(NumericVector x) {
  double ans;
  int n = x.size();
  for (int i = 0; i < n; i++) {
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


/*** R
# pdistR <- function(x, ys)
#  sqrt( (x - ys) ^ 2 )
*/

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
LogicalVector lgl_biggerYC2(NumericVector x, double y) {
  return x > y;
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


// StackOverflow: subset a vector and sort it

// [[Rcpp::export]]
NumericVector rollP3(NumericVector A, int start, int end) {
  NumericVector B = A[seq(start-1, end-1)] ;
  return B.sort() ;
}

// Rcpp gallery: revert a vector

// [[Rcpp::export]]
NumericVector rcppRev(NumericVector x) {
  NumericVector revX = clone<NumericVector>(x);
  std::reverse(revX.begin(), revX.end());
  return revX;
}

////////////////////////////
//     missing vales      //
////////////////////////////

// [[Rcpp::export]]
LogicalVector is_naC(NumericVector x) {
  int n = x.size();
  LogicalVector ans(n);  
  for (int i = 0; i < n; ++i) {
    ans[i] = NumericVector::is_na(x[i]);
  }
  return ans;
}


// [[Rcpp::export]]
LogicalVector is_naC2(NumericVector x) {
  return is_na(x);
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

// [[Rcpp::export]]
NumericVector foov(NumericVector xx, NumericVector yy) {
  NumericVector res1 = ifelse(xx < yy, xx * xx, -(yy * yy));
}

//////////////////////////
//    other examples    //
//////////////////////////

// recursion

// [[Rcpp::export]]
int fC(int n) {
  if (n < 2) return(n);
  else return(fC(n-1) + fC(n-2));
}

// Calling a function

double square( double x ){
  return x*x ;
}

// [[Rcpp::export]]
NumericVector applyC(NumericVector xx){
  return sapply(xx, square);
}

// [[Rcpp::export]]
NumericVector callFunction(NumericVector x,
			   Function f) {
  NumericVector res = f(x);
  return res;
}

// callFunction(x, summary)


// ============================

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

/*** R
library(microbenchmark)
x <- runif(1e5)
microbenchmark(
  mean(x),
  meanC(x)
)
*/


// [[Rcpp::export]]
NumericVector f2(NumericVector x) {
  int n = x.size();
  NumericVector out(n);

  out[0] = x[0];
  for(int i = 1; i < n; ++i) {
    out[i] = out[i - 1] + x[i];
  }
  return out;
}

// [[Rcpp::export]]
bool f3(LogicalVector x) {
  int n = x.size();

  for(int i = 0; i < n; ++i) {
    if (x[i]) return true;
  }
  return false;
}

// [[Rcpp::export]]
int f4(Function pred, List x) {
  int n = x.size();

  for(int i = 0; i < n; ++i) {
    LogicalVector res = pred(x[i]);
    if (res[0]) return i + 1;
  }
  return 0;
}

// [[Rcpp::export]]
NumericVector f5(NumericVector x, NumericVector y) {
  int n = std::max(x.size(), y.size());
  NumericVector x1 = rep_len(x, n);
  NumericVector y1 = rep_len(y, n);

  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    out[i] = std::min(x1[i], y1[i]);
  }

  return out;
}

// [[Rcpp::export]]
NumericVector attribs() {
  NumericVector out = NumericVector::create(1, 2, 3);

  out.names() = CharacterVector::create("a", "b", "c");
  out.attr("my-attr") = "my-value";
  out.attr("class") = "my-class";

  return out;
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List scalar_missings() {
  int int_s = NA_INTEGER;
  String chr_s = NA_STRING;
  bool lgl_s = NA_LOGICAL;
  double num_s = NA_REAL;

  return List::create(int_s, chr_s, lgl_s, num_s);
}


#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double sum4(NumericVector x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

// [[Rcpp::export]]
double sum4x(NumericVector x) {
  return std::accumulate(x.begin(), x.end(), 0);
}
