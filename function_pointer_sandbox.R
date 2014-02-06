setwd("~/Documents/fast_point_mixtures")
library(Rcpp)
library(microbenchmark)
library(ggplot2)



sourceCpp("function_pointer_sandbox.cpp")


N <- 1e2
x <- matrix(1:(2 * N), ncol=2)

foo <- MatrixPointer(x)


microbenchmark(Raw(x))
microbenchmark(AddProduct(x))
microbenchmark(AddMatrixProduct(x))
microbenchmark(AddDoubleProduct(x))



# Holy shit this is slow
Product <- function(x) { x[1] * x[2] }
#microbenchmark(WithFun(x, Product))