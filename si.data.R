library(MASS)
si.data.Wong.Chan <- function(dimension = 10, sample.size = 200) {
  sample.matrix <- matrix(0, sample.size, dimension + 4)
  for (i in 1:sample.size) {
    z <- mvrnorm(1, mu = rep(0, dimension), Sigma = diag(dimension))
    sample.matrix[i, 1] <- exp(z[1] / 2)
    sample.matrix[i, 2] <- z[2] / (exp(z[1]) + 1)
    sample.matrix[i, 3] <- (z[1] * z[3] + 0.6)^(3)
    sample.matrix[i, 4] <- (z[2] + z[4] + 20)^2
    sample.matrix[i, 5:dimension] <- z[5:dimension]
    p <- 1 / (exp(0.5 * z[1] + 0.1 * z[4]) + 1)
    sample.matrix[i, dimension + 1] <- rbinom(1, 1, p)
    temp <- rnorm(1)
    sample.matrix[i, dimension + 2] <- 210 + (1.5 * sample.matrix[i, dimension + 1] - 0.5) * (13.7 * z[1] + 13.7 * z[2] + 13.7 * z[3] + 13.7 * z[4]) + temp
    sample.matrix[i, dimension + 3] <- 210 + 1 * (13.7 * z[1] + 13.7 * z[2] + 13.7 * z[3] + 13.7 * z[4])
    sample.matrix[i, dimension + 4] <- 210 - 0.5 * (13.7 * z[1] + 13.7 * z[2] + 13.7 * z[3] + 13.7 * z[4])
  }
  return(list(covariate = sample.matrix[, 1:dimension], treat = sample.matrix[, dimension + 1], outcome = sample.matrix[, dimension + 2]))
}

si.data.interaction <- function(dimension = 10, sample.size = 200) {
  sample.matrix <- matrix(0, sample.size, dimension + (dimension * (dimension + 1) / 2 + 4))
  covariance.matrix.A <- matrix(0, dimension, dimension)
  covariance.matrix.B <- matrix(0, dimension, dimension)
  for (i in 1:dimension) {
    for (j in 1:dimension) {
      covariance.matrix.A[i, j] <- 2^(-abs(i != j))
    }
  }
  for (i in 1:dimension) {
      covariance.matrix.B[i, i] <- 1
  }
  for (i in 1:sample.size) {
    p <- rbinom(1, 1, 0.5)
    z1 <- mvrnorm(1, mu = rep(1, dimension), Sigma = covariance.matrix.A)
    z0 <- mvrnorm(1, mu = rep(1, dimension), Sigma = covariance.matrix.B)
    if (p == 1) {
      z <- z1
    }
    if (p == 0) {
      z <- z0
    }
    sample.matrix[i, 1:dimension] <- z
    for (f in 1:dimension) {
      for (k in 1:f) {
        sample.matrix[i, dimension + f * (f - 1) / 2 + k] <- z[k] * z[f]
      }
    }
    sample.matrix[i, dimension + (dimension * (dimension + 1) / 2 + 1)] <- p
    temp <- rnorm(1)
    sample.matrix[i, dimension + (dimension * (dimension + 1) / 2 + 3)] <- 1 * sum(z) + 2 * (z[1] * z[2] + z[2] * z[3] + z[3] * z[4] + z[4] * z[5] + z[5] * z[6] + z[6] * z[7] + z[7] * z[8] + z[8] * z[9] + z[9] * z[10] + z[10] * z[1])
    sample.matrix[i, dimension + (dimension * (dimension + 1) / 2 + 4)] <- 1 * sum(z) + (z[1] * z[2] + z[2] * z[3] + z[3] * z[4] + z[4] * z[5] + z[5] * z[6] + z[6] * z[7] + z[7] * z[8] + z[8] * z[9] + z[9] * z[10] + z[10] * z[1])
    sample.matrix[i, dimension + (dimension * (dimension + 1) / 2 + 2)] <- p * sample.matrix[i, (dimension + dimension * (dimension + 1) / 2 + 3)] + (1 - p) * sample.matrix[i, dimension + (dimension * (dimension + 1) / 2 + 4)] + temp
  }
  true_Y <- mean(sample.matrix[i, dimension + (dimension * (dimension + 1) / 2 + 3)] -
                   sample.matrix[i, dimension + (dimension * (dimension + 1) / 2 + 4)])
  return(list(covariate = sample.matrix[, 1:(dimension + dimension * (dimension + 1) / 2)],
              treat = sample.matrix[, dimension + dimension * (dimension + 1) / 2 + 1],
              outcome = sample.matrix[, dimension + dimension * (dimension + 1) / 2 + 2],
              true_Y = true_Y))
}

si.data.extreme.mean.diff <- function(dimension = 10, sample.size = 200) {
  covmatrix <- matrix(0, dimension, dimension)
  treat <- rbinom(sample.size, 1, 0.5)
  z1 <- mvrnorm(sample.size, mu = rep(1, dimension), Sigma = diag(dimension))
  z0 <- mvrnorm(sample.size, mu = rep(0, dimension), Sigma = diag(dimension))
  covariate <- treat * z1 + (1 - treat) * z0
  noise <- rnorm(sample.size)
  outcome <- (1 + treat) * apply(covariate[, 1:10], 1, sum) + noise
  return(list(covariate = covariate, treat = treat, outcome = outcome))
}

si.data.extreme.sample.size <- function(dimension = 10, sample.size = 1000) {
  covariate <- matrix(rnorm(sample.size * dimension) + 1, nrow = sample.size, ncol = dimension)
  ps <- 1 / (1 + 19 * exp(apply(covariate[, 1:10], 1, sum) - 10))
  treat <- rep(NA, sample.size)
  for (i in 1:sample.size) {
    treat[i] <- rbinom(1, 1, ps[i])
  }
  noise <- rnorm(sample.size)
  outcome <- (1 + treat) * apply(covariate[, 1:10], 1, sum) + noise
  return(list(covariate = covariate, treat = treat, outcome = outcome))
}

si.data.heavy.tailed <- function(dimension = 10, sample.size = 200) {
  covariate <- matrix(rnorm(sample.size * dimension), nrow = sample.size, ncol = dimension)
  covariate <- exp(covariate)
  ps <- 1 / (1 + 0.1 * exp(apply(covariate[, 1:5], 1, sum) - 5))
  treat <- rep(NA, sample.size)
  for (i in 1:sample.size) {
    treat[i] <- rbinom(1, 1, ps[i])
  }
  noise <- rnorm(sample.size)
  outcome <- (1 + treat) * apply(covariate[, 1:5], 1, sum) + noise
  true_Y <- mean(apply(covariate[, 1:5], 1, sum))
  return(list(covariate = covariate, treat = treat, outcome = outcome, true_Y = true_Y))
}
