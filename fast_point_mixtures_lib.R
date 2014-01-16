
setwd("~/Documents/fast_point_mixtures")
library(Rcpp)
sourceCpp("fast_point_mixtures.cpp")


REvaluateCDF <- function(x, means, vars, probs) {
  k <- length(means)
  p.values <- sapply(1:k, function(i) {
    return(pnorm(x, mean=means[i], sd=sqrt(vars[i])))
  })
  return(p.values %*% probs)
}


REvaluateDensity <- function(x, means, vars, probs) {
  k <- length(means)
  p.values <- sapply(1:k, function(i) {
    return(dnorm(x, mean=means[i], sd=sqrt(vars[i])))
  })
  return(p.values %*% matrix(probs, ncol=1))
}


EMMixtureOfNormals <- function(x, means, vars,
                               log.prior.probs=rep(-log(length(mean.start))),
                               row.weights=rep(1, length(x)),
                               fit.means=T, fit.vars=F, fit.probs=F,
                               verbose=T, max.iters=1e6,
                               diff.tolerance=1e-4, sd.tolerance=1,
                               num.check.points=min(c(100, length(x))),
                               num.check.buckets=10,
                               check.convergence.every=10) {
  
  if (!any(c(fit.means, fit.vars, fit.probs))) {
    stop("You must fit the means, vars, or probs.")
  }
  #check.points <- quantile(x, 1:num.check.points / (num.check.points + 1))
  check.points <- x
  check.buckets <- c(0, 1:num.check.buckets / num.check.buckets)
  
  # This is the standard deviation of the binomial variance stabilizing transform
  # arcsin(sqrt(p)).
  uniform.bucket.sd <- 1 / sqrt(4 * num.check.points)
  
  # Intialize
  iter <- 0
  k <- length(means)
  probs <- matrix(0, ncol=k, nrow=N)
  cdf <- old.cdf <- REvaluateCDF(check.points, means, true.vars, exp(log.prior.probs))
  converged <- F
  while(!converged) {
    result  <- NormalPointMixtureSummary(x, rep(0, length(x)), means, vars, probs,
                                         log.prior.probs, row.weights)
    if (fit.means) {
      means <- result$x / result$prob_tot
    }
    if (fit.vars) {
      kMinVar <- 1e-6
      vars <- result$x2 / result$prob_tot - (result$x / result$prob_tot)^2
      vars[vars <= 0 | result$prob_tot == 0] <- kMinVar
    }
    if (fit.probs) {
      kMinProb <- 1e-6
      log.prior.probs <- log((result$prob_tot + kMinProb) / sum(result$prob_tot + kMinProb))
    }
    
    # Check for convergence
    if (iter %% check.convergence.every == 0) {
      cdf <- REvaluateCDF(check.points, means, vars, exp(log.prior.probs))
      diff <- mean(abs(cdf - old.cdf))
      old.cdf <- cdf
      if (diff < diff.tolerance) {
        if (verbose) {
          print("Converged due to CDF difference.")
        }
        converged <- T
      } else {
        # Check if the variance of the p values is close to the variance of a
        # uniform distribution.
        cdf.buckets <- cut(cdf, breaks=check.buckets)
        bucket.totals <- tapply(X=cdf.buckets, INDEX=cdf.buckets, FUN=length) / num.check.points
        
        # If there are no elements in a bucket, you get an NA
        bucket.totals[is.na(bucket.totals)] <- 0
        bucket.sd <- sd(asin(sqrt(bucket.totals))) 
        if (bucket.sd < sd.tolerance * uniform.bucket.sd) {
          if (verbose) {          
            print("Converged due to standard deviation.")
          }
          converged <- T
        }
      }
      if (verbose) {
        print(sprintf("%d: Diff: %f, sd ratio: %f", iter, diff, bucket.sd / uniform.bucket.sd))
      }
    }
    iter <- iter + 1
    if (iter > max.iters) {
      if (verbose) {
        print("Converged due to maximum iterations reached.")
      }
      converged <- T
    }
  }
  return(list(probs=probs, means=means, vars=vars, log.prior.probs=log.prior.probs,
              log.prior.probs=log.prior.probs, iter=iter, cdf=cdf))
}