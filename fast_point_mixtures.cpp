#include <Rcpp.h>
#include <assert.h>
using namespace Rcpp;

void RcppAssert(bool expression, const char *message) {
  if (!expression) {
    // I think range_error is not the right thing, but it's not worth
    // fussing about, and I don't really know how c++ exception handling
    // works.
    throw std::range_error(message);
  }
}

// This isn't faster than the R version.
// [[Rcpp::export]]
NumericVector EvaluateCDF(NumericVector x,
                          NumericVector means, NumericVector vars, NumericVector probs) {
                            
  Environment stats("package:stats");
  Function pnorm = stats["pnorm"];

  //Function PNorm("pnorm");
  NumericVector cdf(x.size());
  for (int component_i = 0; component_i < means.size(); component_i++) {
    NumericVector prob = pnorm(x,
                               Named("mean", means(component_i)),
                               Named("sd", sqrt(vars(component_i))));
    cdf = cdf + probs(component_i) * prob;
  }
  return cdf;
}



// [[Rcpp::export]]
List NormalPointMixtureSummary(NumericVector x, NumericVector x_vars,
                               NumericVector means, NumericVector vars,
                               NumericMatrix probs, NumericVector log_prior_probs,
                               NumericVector row_weights) {
  // Returns a matrix of as many rows as the length of <x> and as many columns
  // as the lengths of <means> and <vars>.
  
  int x_length = x.size(); 
  int mean_length = means.size();
  int vars_length = vars.size();
  
  RcppAssert(mean_length == vars_length,
             "The mean and variance must be the same length.");
                          
  RcppAssert(log_prior_probs.size() == mean_length,
             "Wrong length for log_prior_probs.");
  
  RcppAssert(probs.ncol() == mean_length,
             "Wrong number of columns in probs.");

  RcppAssert(probs.nrow() == x_length,
             "Wrong number of rows in probs.");
             
  RcppAssert(row_weights.size() == x_length,
             "Wrong length for row_weights.");
             
  RcppAssert(x_vars.size() == x_length,
             "Wrong length for x_vars.");

  NumericVector row_maxima(x_length);
  NumericVector row_totals(x_length);

  NumericVector x_summary(mean_length);
  NumericVector x2_summary(mean_length);
  NumericVector prob_summary(mean_length);
  
  int row_i = 0;
  int col_i = 0;

  // Calculate the log likelihoods.
  for (row_i = 0; row_i < x_length; row_i++) {
    for (col_i = 0; col_i < mean_length; col_i++) {
      
      // TODO: you can eliminate an entire pass through the matrix by being smarter
      // about the totals.  In this first pass through the columns, exponentiate and
      // keep track of the total.  If you come to a number that will cause an over or
      // underflow, multiply both the total and all future values by a scalar that
      // brings it within a usable range.  Keep track of the scalar and keep modifying
      // it as needed to prevent over and underflows.
  
      // TODO: make this generic with a function pointer.
      probs(row_i, col_i) = log_prior_probs(col_i) -
                            0.5 * pow(x(row_i) - means(col_i), 2) / (vars(col_i) + x_vars(row_i)) -
                            0.5 * log(vars(col_i) + x_vars(row_i));
      if (probs(row_i, col_i) > row_maxima(row_i)) {
        row_maxima(row_i) = probs(row_i, col_i);
      }
    }
  }
  
  // Exponentiate.
  for (row_i = 0; row_i < x_length; row_i++) {
    for (col_i = 0; col_i < mean_length; col_i++) {
      probs(row_i, col_i) = exp(probs(row_i, col_i) - row_maxima(row_i));
      row_totals(row_i) += probs(row_i, col_i);
    }
  }
  
  // Normalize.
  for (row_i = 0; row_i < x_length; row_i++) {
    for (col_i = 0; col_i < mean_length; col_i++) {
      probs(row_i, col_i) = probs(row_i, col_i) / row_totals(row_i);
      x_summary(col_i) +=  row_weights(row_i) * x(row_i) * probs(row_i, col_i);
      x2_summary(col_i) += row_weights(row_i) *
                           (pow(x(row_i), 2) - x_vars(row_i)) * probs(row_i, col_i);
      prob_summary(col_i) += row_weights(row_i) * probs(row_i, col_i);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("x") = x_summary,
                            Rcpp::Named("x2") = x2_summary,
                            Rcpp::Named("prob_tot") = prob_summary);
}


// [[Rcpp::export]]
List PointMixtureSummary(NumericVector x, NumericVector x_vars,
                         NumericVector means, NumericVector vars,
                         NumericMatrix probs, NumericVector log_prior_probs,
                         NumericVector row_weights) {
  // Returns a matrix of as many rows as the length of <x> and as many columns
  // as the lengths of <means> and <vars>.
  
  int x_length = x.size(); 
  int mean_length = means.size();
  int vars_length = vars.size();
  
  RcppAssert(mean_length == vars_length,
             "The mean and variance must be the same length.");
                          
  RcppAssert(log_prior_probs.size() == mean_length,
             "Wrong length for log_prior_probs.");
  
  RcppAssert(probs.ncol() == mean_length,
             "Wrong number of columns in probs.");

  RcppAssert(probs.nrow() == x_length,
             "Wrong number of rows in probs.");
             
  RcppAssert(row_weights.size() == x_length,
             "Wrong length for row_weights.");
             
  RcppAssert(x_vars.size() == x_length,
             "Wrong length for x_vars.");

  NumericVector row_maxima(x_length);
  NumericVector row_totals(x_length);

  NumericVector x_summary(mean_length);
  NumericVector x2_summary(mean_length);
  NumericVector prob_summary(mean_length);
  
  int row_i = 0;
  int col_i = 0;

  // Calculate the log likelihoods.
  for (row_i = 0; row_i < x_length; row_i++) {
    for (col_i = 0; col_i < mean_length; col_i++) {
      
      // TODO: you can eliminate an entire pass through the matrix by being smarter
      // about the totals.  In this first pass through the columns, exponentiate and
      // keep track of the total.  If you come to a number that will cause an over or
      // underflow, multiply both the total and all future values by a scalar that
      // brings it within a usable range.  Keep track of the scalar and keep modifying
      // it as needed to prevent over and underflows.
  
      // TODO: make this generic with a function pointer.
      probs(row_i, col_i) = log_prior_probs(col_i) -
                            0.5 * pow(x(row_i) - means(col_i), 2) / (vars(col_i) + x_vars(row_i)) -
                            0.5 * log(vars(col_i) + x_vars(row_i));
      if (probs(row_i, col_i) > row_maxima(row_i)) {
        row_maxima(row_i) = probs(row_i, col_i);
      }
    }
  }
  
  // Exponentiate.
  for (row_i = 0; row_i < x_length; row_i++) {
    for (col_i = 0; col_i < mean_length; col_i++) {
      probs(row_i, col_i) = exp(probs(row_i, col_i) - row_maxima(row_i));
      row_totals(row_i) += probs(row_i, col_i);
    }
  }
  
  // Normalize.
  for (row_i = 0; row_i < x_length; row_i++) {
    for (col_i = 0; col_i < mean_length; col_i++) {
      probs(row_i, col_i) = probs(row_i, col_i) / row_totals(row_i);
      x_summary(col_i) +=  row_weights(row_i) * x(row_i) * probs(row_i, col_i);
      x2_summary(col_i) += row_weights(row_i) *
                           (pow(x(row_i), 2) - x_vars(row_i)) * probs(row_i, col_i);
      prob_summary(col_i) += row_weights(row_i) * probs(row_i, col_i);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("x") = x_summary,
                            Rcpp::Named("x2") = x2_summary,
                            Rcpp::Named("prob_tot") = prob_summary);
}
