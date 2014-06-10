library(RUnit)
kRepoLocation <- Sys.getenv("GIT_REPO_LOC")

source(file.path(kRepoLocation, "fast_point_mixtures/fast_point_mixtures_lib.R"))


RawLikelihoodLogLik <- function(x, x.vars, means, vars) {
  n <- length(x)
  k <- length(means)
  StopIfNotEqual(n, length(x.vars), "x.vars is the wrong length")  
  StopIfNotEqual(k, length(vars), "vars is the wrong length")
    RowLogLik <- function(i) {
    return(dnorm(x=x[i], mean=means, sd=(sqrt(vars + x.vars[i])), log=TRUE))
  }
  return(do.call(rbind, lapply(1:n, RowLogLik)))
}


testNormalPointLogLik <- function() {
  k <- 2
  k.seq <- (1:k) / k
  means <- k * k.seq
  vars <- (0.7^2) * k.seq
  
  x <- rep(means, each=2) + 0.1 * (1:n) / n
  n <- length(x)
  x.vars <- rep(0.1 ^ 2, n)
  
  prior.probs <- (k.seq + 2)
  prior.probs <- prior.probs / sum(prior.probs)
  log.prior.probs <- log(prior.probs)
  
  probs <- matrix(0, nrow=n, ncol=k)
  row.maxima <- NormalLogLikelihoodMatrix(x, x.vars,
                                          means, vars,
                                          probs, log.prior.probs)

  raw.log.lik <- RawLikelihoodLogLik(x=x, x.vars=x.vars, means=means, vars=vars)
  raw.log.probs <- raw.log.lik + rep(log.prior.probs, each=n)
  
  # The two may be off by a constant.
  differences <- as.numeric(raw.log.probs - probs)
  checkEqualsNumeric(0, abs(min(differences) - max(differences)),
                     "log likelihood matrix incorrect")
  
  checkEqualsNumeric(apply(probs, 1, max), row.maxima,
                     "Row maxima are incorrect.")
}


testNormalPointMixtureFun <- function() {
  k <- 2
  k.seq <- (1:k) / k
  means <- k * k.seq
  vars <- (0.6^2) * k.seq
  
  x <- rep(means, each=2) + 0.1 * (1:n) / n
  n <- length(x)
  x.vars <- rep(0.1 ^ 2, n)
  
  prior.probs <- (k.seq + 2)
  prior.probs <- prior.probs / sum(prior.probs)
  log.prior.probs <- log(prior.probs)
  row.weights <- ((1:n) / n + 10)
  row.weights <- row.weights / sum(row.weights)
  
  probs <- matrix(0, nrow=n, ncol=k)
  probs.summary <- NormalPointMixtureSummary(x, x.vars,
                                             means, vars,
                                             probs, log.prior.probs,
                                             row.weights)

  raw.log.lik <- RawLikelihoodLogLik(x=x, x.vars=x.vars, means=means, vars=vars)
  raw.log.probs <- raw.log.lik + rep(log.prior.probs, each=n)
  raw.log.probs <- raw.log.probs - apply(raw.log.probs, 1, max)
  raw.probs <- exp(raw.log.probs)
  raw.probs <- raw.probs / apply(raw.probs, 1, sum)
  
  checkEqualsNumeric(raw.probs, probs, "Probs incorrect")
  checkEqualsNumeric(colSums(row.weights * raw.probs), probs.summary$prob_tot,
                     "prob_tot incorrect")
  checkEqualsNumeric(colSums(row.weights * x * raw.probs), probs.summary$x,
                     "x summary incorrect")
  checkEqualsNumeric(colSums(row.weights * (x ^ 2 - x.vars) * raw.probs), probs.summary$x2,
                     "x^2 summary incorrect")
}
