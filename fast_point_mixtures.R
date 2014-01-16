setwd("~/Documents/fast_point_mixtures")
library(Rcpp)
library(microbenchmark)
library(ggplot2)
sourceCpp("fast_point_mixtures.cpp")

N <- 10000
k <- 10
true.means <- 1:k
true.vars <- rep(0.1^2, k)
true.probs <- rep(1, k)
true.probs <- true.probs / sum(true.probs)
#noise.var <- rep(0.2^2, N)
noise.var <- rep(0, N)

components <- t(rmultinom(N, prob=true.probs, size=1))
x.means <- components %*% matrix(true.means)
x.vars <- components %*% matrix(true.vars) + noise.var
x <- rnorm(N, mean=x.means, sd=sqrt(x.vars))
var(x[x.means == 1])

#ggplot() +
#  geom_line(aes(x=x, y=REvaluateDensity(x, true.means, true.vars, true.probs), color="true"), lwd=2) +
#  geom_density(aes(x=x, color="empirical"), adjust=0.2)

row.weights <- rep(1, N)
log.prior.probs <- log(true.probs)

# Starting value for the means
means <- quantile(x, (1:k) / (k + 1))
#means <- 2 * 1:10
results <- EMMixtureOfNormals(x=x, means=true.means, vars=rep(10, k),
                              fit.mean=F, fit.vars=T, fit.probs=F,
                              max.iters=100,
                              check.convergence.every=10,
                              num.check.buckets=N / 100,
                              num.check.points=N,
                              log.prior.probs=log.prior.probs)

print(results$vars)
qplot(results$cdf, geom="histogram")
ggplot() +
  geom_line(aes(x=x, y=REvaluateDensity(x, true.means, true.vars, true.probs), color="true"), lwd=2) +
  geom_line(aes(x=x, y=REvaluateDensity(x, results$means, results$vars, true.probs), color="fit"), lwd=2) +
  ylim(0, 1)


true.p.values <- REvaluateCDF(x, means=true.means, vars=true.vars, true.probs)
p.values <- REvaluateCDF(x, results$means, results$vars, true.probs)
ggplot() +
  geom_density(aes(x=p.values, color="fit"), adjust=0.2) +
  geom_density(aes(x=true.p.values, color="true"), adjust=0.2)


hist(p.values)
ggplot() +
  geom_density(aes(x, color="all")) +
  geom_density(aes(x[p.values < 0.01], color="low.p")) +
  geom_density(aes(x[p.values > 1 - 0.01], color="high.p"))  


####3

x <- rnorm(N, mean=1, sd=2)
results <- EMMixtureOfNormals(x=x, means=1, vars=0.1,
                              fit.mean=F, fit.vars=T, fit.probs=F,
                              max.iters=200,
                              sd.tolerance=0.01,
                              check.convergence.every=1,
                              num.check.points=N,
                              log.prior.probs=0)

