#include <Rcpp.h>
using namespace Rcpp;

/* Doesn't work
// [[Rcpp::export]]
NumericVector RawDouble(NumericVector x_mat, int x_length) {
  double x[x_length] = x_mat;
  NumericVector result(1);
  for (int i = 0; i < x_length; i++) {
    result = result + x[i];
  }
  return result;
}
*/


NumericMatrix MatrixPointerFirstElement(NumericMatrix* x_pointer) {
  NumericMatrix x = *x_pointer;
  return (x)(1,1);
}

// [[Rcpp::export]]
NumericMatrix MatrixPointer(NumericMatrix x) {
  return MatrixPointerFirstElement(&x);
}





// [[Rcpp::export]]
NumericVector Raw(NumericMatrix x) {
  NumericVector result(1);
  for (int i = 0; i < x.rows(); i++) {
    // It's the allocation of a new numeric vector, not the function
    // pointer, that slows down the other one.
    //NumericVector y = x(i, _);
    result = result + (x(i, 0) * x(i, 1));
  }
  return result;
}

//////////////////////////////
double CProduct(NumericVector x) {
  return x(0) * x(1);
}

NumericVector WithFun(NumericMatrix x, double (*CombineFun)(NumericVector z)) {
  NumericVector result(1);
  for (int i = 0; i < x.rows(); i++) {
    result = result + CombineFun(x(i, _));
  }
  return result;
}

// This is about one fifth the speed.
// [[Rcpp::export]]
NumericVector AddProduct(NumericMatrix x) {
  return WithFun(x, CProduct);
}


//////////////////////////////
double CMatrixProduct(NumericMatrix x, int i) {
  return x(i, 0) * x(i, 1);
}

NumericVector WithMatrixFun(NumericMatrix x, double (*CombineFun)(NumericMatrix, int)) {
  NumericVector result(1);
  for (int i = 0; i < x.rows(); i++) {
    result = result + CombineFun(x, i);
  }
  return result;
}

// This is about half the speed.
// [[Rcpp::export]]
NumericVector AddMatrixProduct(NumericMatrix x) {
  return WithMatrixFun(x, CMatrixProduct);
}



//////////////////////////////
double CDoubleProduct(double a, double b) {
  return a * b;
}

NumericVector WithDoubleFun(NumericMatrix x, double (*CombineFun)(double, double)) {
  NumericVector result(1);
  for (int i = 0; i < x.rows(); i++) {
    result = result + CombineFun(x(i,0), x(i,1));
  }
  return result;
}

// This is about half the speed.
// [[Rcpp::export]]
NumericVector AddDoubleProduct(NumericMatrix x) {
  return WithDoubleFun(x, CDoubleProduct);
}




/*
//////////////////////////////
// This is crazy slow
// [[Rcpp::export]]
NumericVector RWithFun(NumericMatrix x, Function Combine) {
  NumericVector result(1);
  for (int i = 0; i < x.rows(); i++) {
    NumericVector x_combine = Combine(x(i, _));
    result = result + x_combine;
  }
  return result;
}
*/




/* This doesn't make it faster
//////////////////////////////
double SEXPCProduct(SEXP x_sexp, int i) {
  NumericMatrix x = x_sexp;
  return x(i, 0) * x(i, 1);
}

NumericVector SEXPWithFun(SEXP x_sexp, double (*CombineFun)(NumericMatrix, int)) {
  NumericVector result(1);
  NumericMatrix x = x_sexp;
  for (int i = 0; i < x.rows(); i++) {
    result = result + CombineFun(x_sexp, i);
  }
  return result;
}

// This is about half the speed.
// [[Rcpp::export]]
NumericVector SEXPAddProduct(NumericMatrix x) {
  return WithFun(x, CProduct);
}
*/