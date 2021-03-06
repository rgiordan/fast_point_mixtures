setwd("~/Documents/fast_point_mixtures")
library(Rcpp)
library(microbenchmark)
library(ggplot2)
library(gtools) # for rdirichlet

sourceCpp("fast_point_mixtures.cpp")

N <- 1e2
k <- 3
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


means <- true.means
vars <- true.vars
log.prior.probs <- log(true.probs)
probs <- matrix(0, N, k)

microbenchmark(results <- NormalPointMixtureSummary(x=x, x_vars=rep(0, N),
                                     means=means, vars=vars, probs=probs,
                                     log_prior_probs=log.prior.probs,
                                     row_weights=rep(1, N)))

observations <- cbind(x, x^2)
parameters <- cbind(means, vars)
microbenchmark(results.fp <- NormalPointMixtureSummaryFP(observations=observations,
                                                         parameters=parameters,
                                                         probs=probs,
                                                         log_prior_probs=log.prior.probs,
                                                         row_weights=rep(1, N)))


#############
# Look at the normal version


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


################################
# Dirichlet

N <- 1e4
k <- 3
components <- 2
alpha <- matrix(c(1, 2, 3,
                  5, 4, 3), nrow=components)
alpha.mat <- matrix(alpha, nrow=components)
probs <- (1:components) / sum(1:components)
x.counts <- rmultinom(1, size=N, prob=probs)
x <- do.call(rbind, lapply(1:components, function(i) { rdirichlet(n=x.counts[i],
                                                                  alpha=alpha[i, ])}))
x.components <- do.call(c, lapply(1:components, function(i) { rep(i, x.counts[i]) }))



# Check the RCPP function
i <- 2
log(ddirichlet(x[i,], alpha[1,,drop=F]))
DirichletLogLikelihood(x[i,], alpha[1,])

# Check the fitting function
FitDirichletFromSufficientStats(colSums(log(x)), N, rep(1, k))
rowSums(t(alpha) %*% probs)

# Check the mixture code
probs <- matrix(0, nrow=N, ncol=components)
x.augmented <- cbind(x, log(x))
result  <- DirichletPointMixtureSummaryFP(observations=x.augmented,
                                          alpha=alpha.mat,
                                          probs=probs,
                                          log_prior_probs=rep(0, k),
                                          row_weights=rep(1, N))
result$summary_stats / result$prob_tot



