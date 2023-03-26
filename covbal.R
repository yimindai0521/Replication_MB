# Implementation of covariate balancing method
library(WeightIt)
library(MASS)
library(ACBalancing)
library(ATE.ncb)
library(sbw)

covbal <- function(covariate, treat, outcome, kernel = FALSE, para.method = FALSE) {
  covariate <- scale(as.matrix(covariate))
  treat <- as.vector(treat)
  outcome <- as.vector(outcome)
  data.matrix <- data.frame(covariate, factor(treat))
  sample.size <- dim(covariate)[1]
  dimension <- dim(covariate)[2]
  name.method <- c("unad", "ps", "ebal", "ebcw", "cbps", "MB2", "MB", "sbw", "kernel")
  for (j in 1:length(name.method)) {
    assign(paste(name.method[j], ".", "weight", sep = ""), rep(NA, sample.size))
  }
  unad.weight <- rep(1, sample.size)

  if (para.method == TRUE) {
    character <- names(data.matrix)
    for (j in 1:(dimension + 1)) {
      character[j] <- paste(character[j])
    }
    myformula <- as.formula(paste(character[1 + dimension], paste(" ~ ", paste(character[1:dimension], collapse = "+")))) # initialize formula for weightit function

    # Implementation of different weighting methods
    tryCatch(ps.weight <- weightit(myformula, data = data.matrix, estimand = "ATE", method = "ps")$weights,
      error = function(e) {
        skip_to_next <<- FALSE
      }
    ) # Inverse probability weighting
    ebal.weight <- weightit(myformula, data = data.matrix, estimand = "ATE", method = "ebal")$weights # entropy balancing
    tryCatch(ebcw.weight <- weightit(myformula, data = data.matrix, estimand = "ATE", method = "ebcw")$weights,
      error = function(e) {
        skip_to_next <<- FALSE
      }
    ) # Empirical balancing calibration weighting
    cbps.weight <- weightit(myformula, data = data.matrix, estimand = "ATE", method = "cbps", over = FALSE)$weights # covariate balancing propensity score
    # energy.weight <- weightit(myformula, data = data.matrix, estimand = "ATE", method = "energy")$weights # energy balancing

    # univariate balancing
    bal <- list()
    bal$bal_gri <- c(1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
    bal$bal_cov <- character[1:dimension]
    sbw.result <- FALSE
    while (sum(dim(as.matrix(sbw.result))) != 15) {
      sbw.result <- tryCatch(sbw.result <- sbw(
        dat = data.matrix,
        ind = character[dimension + 1], bal = bal, out = outcome, par = list(par_est = "ate")
      ), error = function(e) {
        skip_to_next <<- FALSE
      })
      bal$bal_gri <- bal$bal_gri[-1]
    }
    if (sum(dim(as.matrix(sbw.result))) == 15) {
      sbw.weight <- sbw.result$dat_weights[, dimension + 2]
    }

    # Mahalanobis balancing
    MB.weight <- rep(NA, sample.size)
    MB1.result <- MB(covariate = covariate, treat = treat, group1 = 0,
                     outcome = rep(0, sample.size), delta.space = c(1e-1, 1e-2, 1e-3, 1e-4),
                     method = "MB", opti.method = "BFGS", rate = 10)
    MB2.result <- MB(covariate = covariate, treat = treat, group1 = 1,
                     outcome = rep(0, sample.size), method = "MB", delta.space = c(1e-1, 1e-2, 1e-3, 1e-4),
                     opti.method = "BFGS", rate = 10)
    MB.weight[treat == 0] <- MB1.result$weight
    MB.weight[treat == 1] <- MB2.result$weight

    MB2.weight <- rep(NA, sample.size)
    MB1.result <- MB(covariate = covariate, treat = treat, group1 = 0,
                     outcome = rep(0, sample.size), delta.space = c(1e-1, 1e-2, 1e-3, 1e-4),
                     method = "MB2", opti.method = "BFGS", rate = 10)
    MB2.result <- MB(covariate = covariate, treat = treat, group1 = 1,
                     outcome = rep(0, sample.size), method = "MB2", delta.space = c(1e-1, 1e-2, 1e-3, 1e-4),
                     opti.method = "BFGS", rate = 10)
    MB2.weight[treat == 0] <- MB1.result$weight
    MB2.weight[treat == 1] <- MB2.result$weight

    # Create weight space
    if (!exists("weight.space")) {
      weight.space <- rbind(
        unad = unad.weight,
        ps = ps.weight, ebal = ebal.weight, cbps = cbps.weight, ebcw = ebcw.weight,
        MB = MB.weight, MB2 = MB2.weight, sbw = sbw.weight
      )
    }
  }

  # kernel-based covariate balancing
  if (kernel == TRUE) {
    Xstd <- transform.sob(covariate)$Xstd # standardize covariate to [0,1]^p
    K <- getGram(Xstd) # get Gram matrix using Sobolev kernel
    nlam <- 50 # design a grid for the tuning parameter
    lams <- exp(seq(log(1e-8), log(1), len = nlam)) # design a grid for the tuning parameter
    fit1 <- ATE.ncb.SN(treat, K, lam1s = lams) # compute weights for T=1
    if (sum(fit1$warns)) cat("lambda bound warning!\n")
    fit0 <- ATE.ncb.SN(1 - treat, K, lam1s = lams) # compute weights for T=0
    if (sum(fit0$warns)) cat("lambda bound warning!\n")

    #### ATE ####
    kernel.ate <- mean(fit1$w * outcome - fit0$w * outcome) # ATE estimate based on kernel-based estimation (truth=10)
    kernel.weight <- (fit1$w + fit0$w) / sample.size

    # Create weight space
    if (!exists("weight.space")) {
      weight.space <- as.matrix(rbind(unad = unad.weight, kernel = kernel.weight))
    } else {
      weight.space <- as.matrix(rbind(weight.space, kernel = kernel.weight))
    }
  }

  # Scaling weight for each weight
  weight.space[, treat == 1] <- weight.space[, treat == 1] / rowSums(weight.space[, treat == 1])
  weight.space[, treat == 0] <- weight.space[, treat == 0] / rowSums(weight.space[, treat == 0])

  # Calculating balancing measure and causal effect of interest
  if (!is.null(weight.space)) {
    ate <- crossprod(t(weight.space[, treat == 1]), outcome[treat == 1]) -
      crossprod(t(weight.space[, treat == 0]), outcome[treat == 0])
    TASMD <- abs(crossprod(t(weight.space[, treat == 1]), covariate[treat == 1, ])) +
      abs(crossprod(t(weight.space[, treat == 0]), covariate[treat == 0, ]))
    GMIM <- rowSums((crossprod(t(weight.space[, treat == 1]), covariate[treat == 1, ]))^2) +
      rowSums((crossprod(t(weight.space[, treat == 0]), covariate[treat == 0, ]))^2)
  }

  # Return value
  return(list(
    weight = weight.space,
    ate = ate,
    tasmd = TASMD,
    GMIM = GMIM
  ))
}
