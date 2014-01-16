
N <- 1000
p <- 0.99
sims <- 1000
x <- rep(NA, sims)
for (i in 1:sims) {
  x[i] <- sum(runif(N) < p) / N
}

log(p * (1 - p) / N)
log(var(x))

log(1 / (4 * N))
log(var(asin(sqrt(x))))