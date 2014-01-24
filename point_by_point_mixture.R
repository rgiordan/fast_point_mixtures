setwd("~/Documents/fast_point_mixtures")
library(Rcpp)
library(microbenchmark)
library(ggplot2)
sourceCpp("fast_point_mixtures.cpp")
source("fast_point_mixtures_lib.R")

N <- 10000
true.k <- 10
true.means <- 1:true.k
true.vars <- rep(0.5^2, true.k)
true.probs <- rep(1, true.k)
true.probs <- true.probs / sum(true.probs)

components <- t(rmultinom(N, prob=true.probs, size=1))
x.means <- components %*% matrix(true.means)

#noise.var <- rep(0.2^2, N)
#noise.var <- rep(0, N)
#x.vars <- components %*% matrix(true.vars) + noise.var

x.vars <- components %*% matrix(true.vars)
x <- rnorm(N, mean=x.means, sd=sqrt(x.vars))
var(x[x.means == 1])
hist(x)

num.check.buckets <- N / 100
check.buckets <- c(0, 1:num.check.buckets / num.check.buckets)
uniform.bucket.sd <- 1 / sqrt(4 * N)

true.cdf <- REvaluateCDF(x=x, means=true.means, vars=true.vars, probs=true.probs)
true.cdf.buckets <- cut(true.cdf, breaks=check.buckets)
true.bucket.totals <- tapply(X=true.cdf.buckets, INDEX=true.cdf.buckets, FUN=length) / N
#plot(true.bucket.totals)
hist(x, 1000)


x <- exp(rnorm(N, mean=2, sd=1))
hist(x)
means <- quantile(x, 1:true.k / (true.k + 1))
vars <- rep(var(x) / true.k, true.k) 
log.prior.probs <- -log(1:true.k)
results <- EMMixtureOfNormals(x=x, means=means, vars=vars,
                              fit.mean=T, fit.vars=T, fit.probs=T,
                              max.iters=1e4,
                              check.convergence.every=10,
                              num.check.buckets=N / 100,
                              num.check.points=N)

cdf.buckets <- cut(results$cdf, breaks=check.buckets)
bucket.totals <- tapply(X=cdf.buckets, INDEX=cdf.buckets, FUN=length) / N

# If there are no elements in a bucket, you get an NA
bucket.totals[is.na(bucket.totals)] <- 0
bucket.sd <- sd(asin(sqrt(bucket.totals)))
deviations <- asin(sqrt(bucket.totals)) - 1 / num.check.buckets
hist(results$cdf)

points <- seq(min(x), max(x), length.out=1e3)
fit.distribution <- REvaluateDensity(x=points, means=results$means, vars=results$vars,
                                     probs=exp(results$log.prior.probs))
ggplot() + 
  geom_line(aes(x=points, y=fit.distribution, color="Fit"), lwd=2) +
  geom_density(aes(x=x, color="Data"), lwd=2)

#########
# Initialize
k <- 1
means <- mean(x)
vars <- var(x)
log.prior.probs <- log(1)

converged <- F
iters <- 10

while (!converged) { 
  results <- EMMixtureOfNormals(x=x, means=means, vars=vars,
                                fit.mean=F, fit.vars=F, fit.probs=T,
                                max.iters=1e4,
                                check.convergence.every=10,
                                num.check.buckets=N / 100,
                                num.check.points=N)
  log.prior.probs <- results$log.prior.probs
  if (exp(log.prior.probs[k]) < 0.0001) { 
    converged <- T
    break
  }
  
  hist(results$cdf)
  cdf.buckets <- cut(results$cdf, breaks=check.buckets)
  bucket.totals <- tapply(X=cdf.buckets, INDEX=cdf.buckets, FUN=length) / N
  
  # If there are no elements in a bucket, you get an NA
  bucket.totals[is.na(bucket.totals)] <- 0
  bucket.sd <- sd(asin(sqrt(bucket.totals)))
  deviations <- asin(sqrt(bucket.totals)) - 1 / num.check.buckets
  ggplot() +
    geom_point(aes(x=check.buckets[-1], y=bucket.totals)) +
    geom_line(aes(x=check.buckets[-1], y=deviations))
  
  if (all(abs(deviations) < bucket.sd | k == iters)) {
    converged <- T
  } else {
    
    #ggplot() + geom_point(aes(x=levels(cdf.buckets), y=bucket.totals)) + geom_hline(aes(x=1 / num.check.buckets))
    problem.x <- x[cdf.buckets == names(which.max(deviations))]
    
    means <- c(results$means, mean(problem.x))
    vars <- c(results$vars, var(problem.x))
    k <- k + 1
    
    print(rbind(means, vars, exp(log.prior.probs)))
    print("---------")
  }
}

