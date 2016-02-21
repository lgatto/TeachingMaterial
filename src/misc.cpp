#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
int test(int x) { 
  return x*x;
}







// [[Rcpp::export]]
bool any_naC(NumericVector x) {
  return is_true(any(is_na(x)));
}

// [[Rcpp::export]]
bool any_naC2(NumericVector x) {
  return any(is_na(x));
}

// [[Rcpp::export]]
LogicalVector any_naC3(NumericVector x) {
  return any(is_na(x));
}


// [[Rcpp::export]]
List missing_sampler() {
  return(List::create(
		      NumericVector::create(NA_REAL), 
		      IntegerVector::create(NA_INTEGER),
		      LogicalVector::create(NA_LOGICAL), 
		      CharacterVector::create(NA_STRING)));
}

// [[Rcpp::export]]
List testlist() {
  return(List::create(
		      NumericVector::create(1.0, 10.0), 
		      IntegerVector::create(NA_INTEGER),
		      LogicalVector::create(NA_LOGICAL), 
		      CharacterVector::create(NA_STRING)));
}

// [[Rcpp::export]]
List testlist2() {
  NumericVector x(10);
  return(List::create(
		      x, 
		      IntegerVector::create(NA_INTEGER),
		      LogicalVector::create(NA_LOGICAL), 
		      CharacterVector::create(NA_STRING)));
}

// [[Rcpp::export]]
IntegerVector abc() {
  IntegerVector x = seq_len(10);  
  return x[seq(1, 4)];
}

