




.f <- function(x, use.log=T) {
  log.f <- -x - 2 * log(x)
  if (!use.log) {
    return(ifelse(x <= 0, 0, exp(log.f)))
  } else {
    return(ifelse(x <= 0, 0, log.f))
  }
}


f <- function(x, use.log=T) {
  log.f <- dgamma(x, shape=shape, rate=rate, log=T)
  if (!use.log) {
    return(ifelse(x <= 0, 0, exp(log.f)))
  } else {
    return(ifelse(x <= 0, -Inf, log.f))
  }
}

norm.f <- function(x, m, use.log=T) {
  return(dnorm(x, mean=m, sd=sigma, log=T))
}


#x <- seq(0, 3, length.out=1e4)
#plot(x, f(x))

shape <- 3
rate <- 3
k <- 10
m.start <- seq(1, 3, length.out=k)

m <- m.start
xt <- m.start

sigma <- 0.001

for (i in 1:5000) {  
  # Each column is a component, each row a point.
  log.probs <- do.call(cbind, lapply(m, function(z) { norm.f(xt, z) }))
  probs <- rowSums(exp(log.probs - max(log.probs)))
  probs <- probs / sum(probs)
  log.weights <- f(xt) - log(probs)
  weights <- exp(log.weights - max(log.weights))
  weights <- weights / sum(weights)
  
  # Use monte carlo to calculate the quantiles
  sims <- 1e4
  components <- rmultinom(1, sims, weights)
  x.sim <- unlist(lapply(1:length(components), function(i) { rnorm(n=components[i],
                                                                   mean=m[i],
                                                                   sd=sigma) }))
  print(sprintf("%f, %f", mean(x.sim), var(x.sim)))
  m <- quantile(x.sim, (1:k) / (k + 1))
  names(m) <- NULL
  xt <- m
  #print(weights)
}

true.sim <- rgamma(n=sims, rate=rate, shape=shape)
mean(true.sim)
var(true.sim)

library(ggplot2)
ggplot() +
  geom_density(aes(x=true.sim, color="true"), lwd=3) +
  geom_density(aes(x=x.sim, color="fit"), lwd=3)

quantile(x.sim, (1:k) / (k + 1))
quantile(true.sim, (1:k) / (k + 1))
