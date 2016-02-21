library("Rcpp")

add_one <- cfunction(c(x = "numeric"), '
  REAL(x)[0] = REAL(x)[0] + 3;
  return(x);
')
x <- 1
y <- x
add_one(x)
x
y

add_two <- cfunction(c(x = "numeric"), '
  SEXP x_copy;
  PROTECT(x_copy = duplicate(x));
  REAL(x_copy)[0] = REAL(x_copy)[0] + 4;
  UNPROTECT(1);
  return(x_copy);
')
x <- 1
y <- x
add_two(x)
x
y


## see NumericVector zs = clone(ys) in Rcpp



cppFunction('
  int add(int x, int y, int z) {
    int sum = x + y + z;
    return sum;
  }'
)

add # like a regular R function, printing displays info about the function

add(1, 2, 3)


signR <- function(x) {
  if (x > 0) {
    1
  } else if (x == 0) {
    0
  } else {
    -1
  }
}

cppFunction('
int signC(int x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}
')



sumR <- function(x) {
  total <- 0
  for (i in seq_along(x)) {
    total <- total + x[i]
  }
  total
}

cppFunction('
double sumC(NumericVector x) {
  double total = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    total += x[i];
  }
  return total;
}
')


library(microbenchmark)
x <- runif(1e3)
microbenchmark(
    sum(x),
    sumR(x),
    sumC(x)
  )



pdistR <- function(x, ys) {
    sqrt( (x - ys) ^ 2 )
  }


cppFunction('
NumericVector pdistC(double x, DoubleVector ys) {
  int n = x.size();
  NumericVector out(n);
  // or NumericVector out = clone(ys);
  for (int i = 0; i < n; i++) {
    out[i] = sqrt( pow(x - ys[i], 2) );
  }
  return out;
}
')

## See sugar for Rcpp vectorisation
