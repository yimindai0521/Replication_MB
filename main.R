# main function
source("si.data.R")
source("covbal.R")

iteration <- 1000
sample.size.array <- c(200, 200, 200, 1000, 200)
dimension.array <- c(10, 65, 10, 10, 10)
name.method <- c("unad", "ps", "ebal", "cbps", "ebcw", "MB", "MB2", "sbw", "kernel")

for (i in 1:5) {
  assign(paste("ate", i, sep = ""), matrix(0, iteration, length(name.method)))
  assign(paste("weight", i, sep = ""), array(0, dim = c(iteration, length(name.method), sample.size.array[i])))
  assign(paste("datasummary", i, sep = ""), array(0, dim = c(iteration, dimension.array[i] + 2, sample.size.array[i])))
  assign(paste("TASMDsummary", i, sep = ""), array(0, dim = c(iteration, dimension.array[i], length(name.method))))
  assign(paste("GMIM", i, sep = ""), matrix(0, iteration, length(name.method)))
}

data2Y <- rep(NA, iteration)
data5Y <- rep(NA, iteration)

set.seed(1124)
for (i in 1:iteration) {
  # Generating data
  data1 <- si.data.Wong.Chan()
  result <- covbal(data1$covariate, data1$treat, data1$outcome, kernel = TRUE, para.method = TRUE)
  ate1[i, ] <- result$ate
  weight1[i, , ] <- result$weight
  datasummary1[i, , ] <- cbind(data1$covariate, data1$treat, data1$outcome)
  TASMDsummary1[i, , ] <- t(result$tasmd)
  GMIM1[i, ] <- result$GMIM
}
ate1.result <- colMeans(ate1)
names(ate1.result) <- name.method
round(ate1.result - 0, 2)
sd.ate1 <- apply(ate1, 2, sd)
round(sd.ate1, 2)
rmse1 <- (colMeans((ate1 - 0)^2))^{1/2}
round(rmse1, 2)
round(colMeans(colMeans(TASMDsummary1)), 2)
GMIM1.result <- colMeans(GMIM1)
round(GMIM1.result, 2)



set.seed(1124)
for (i in 1:iteration) {
  data2 <- si.data.interaction()
  result1 <- covbal(data2$covariate, data2$treat, data2$outcome, kernel = FALSE, para.method = TRUE)
  result2 <- covbal(data2$covariate[, 1:10], data2$treat, data2$outcome, kernel = TRUE, para.method = FALSE)
  ate2[i, ] <- c(result1$ate, result2$ate[2])
  weight2[i, , ] <- rbind(result1$weight, result2$weight[2, ])
  datasummary2[i, , ] <- cbind(data2$covariate, data2$treat, data2$outcome)
  kernel.tasmd <- as.vector(abs(crossprod(result2$weight[2, data2$treat == 1], scale(data2$covariate)[data2$treat == 1, ])) +
    abs(crossprod(result2$weight[2, data2$treat == 0], scale(data2$covariate)[data2$treat== 0, ])))
  TASMDsummary2[i, , ] <- t(rbind(result1$tasmd, kernel.tasmd))
  kernel.GMIM <- rowSums((crossprod(t(result2$weight[, data2$treat == 1]), data2$covariate[data2$treat == 1, ]))^2) +
    rowSums((crossprod(t(result2$weight[, data2$treat == 0]), data2$covariate[data2$treat == 0, ]))^2)
  GMIM2[i, ] <- c(result1$GMIM, kernel.GMIM[2])
  data2Y[i] <- data2$true_Y
}
ate2.result <- colMeans(ate2)
names(ate2.result) <- name.method
round(ate2.result - mean(data2Y), 2)
sd.ate2 <- apply(ate2, 2, sd)
round(sd.ate2, 2)
rmse2 <- (colMeans((ate2 - mean(data2Y))^2))^{1/2}
round(rmse2, 2)
round(colMeans(colMeans(TASMDsummary2)), 2)
GMIM2.result <- colMeans(GMIM2)
round(GMIM2.result, 2)

set.seed(1124)
for (i in 1:iteration) {
  data3 <- si.data.extreme.mean.diff()
  result <- covbal(data3$covariate, data3$treat, data3$outcome, kernel = TRUE, para.method = TRUE)
  ate3[i, ] <- result$ate
  weight3[i, , ] <- result$weight
  datasummary3[i, , ] <- cbind(data3$covariate, data3$treat, data3$outcome)
  TASMDsummary3[i, , ] <- result$tasmd
  GMIM3[i, ] <- result$GMIM
}
ate3.result <- colMeans(ate3)
names(ate3.result) <- name.method
round(ate3.result - 5, 2)
sd.ate3 <- apply(ate3, 2, sd)
round(sd.ate3, 2)
rmse3 <- (colMeans((ate3 - 5)^2))^{1/2}
round(rmse3, 2)
round(colMeans(colMeans(TASMDsummary3)), 2)
GMIM3.result <- colMeans(GMIM3)
round(GMIM3.result, 2)

set.seed(1124)
for (i in 1:iteration) {
  data4 <- si.data.extreme.sample.size()
  result <- covbal(data4$covariate, data4$treat, data4$outcome, kernel = TRUE, para.method = TRUE)
  ate4[i, ] <- result$ate
  weight4[i, , ] <- result$weight
  datasummary4[i, , ] <- cbind(data4$covariate, data4$treat, data4$outcome)
  TASMDsummary4[i, , ] <- t(result$tasmd)
  GMIM4[i, ] <- result$GMIM
}
ate4.result <- colMeans(ate4)
names(ate4.result) <- name.method
round(ate4.result - 10, 2)
sd.ate4 <- apply(ate4, 2, sd)
round(sd.ate4, 2)
rmse4 <- (colMeans((ate4 - 10)^2))^{1/2}
round(rmse4, 2)
round(colMeans(colMeans(TASMDsummary4)), 2)
GMIM4.result <- colMeans(GMIM4)
round(GMIM4.result, 2)

set.seed(1124)
for (i in 1:iteration) {
  data5 <- si.data.heavy.tailed()
  result <- covbal(data5$covariate, data5$treat, data5$outcome, kernel = TRUE, para.method = TRUE)
  ate5[i, ] <- result$ate
  weight5[i, , ] <- result$weight
  datasummary5[i, , ] <- cbind(data5$covariate, data5$treat, data5$outcome)
  TASMDsummary5[i, , ] <- t(result$tasmd)
  GMIM5[i, ] <- result$GMIM
  data5Y[i] <- data5$true_Y
}
ate5.result <- colMeans(ate5)
names(ate5.result) <- name.method
round(ate5.result - mean(data5Y), 2)
sd.ate5 <- apply(ate5, 2, sd)
round(sd.ate5, 2)
rmse5 <- (colMeans((ate5 - mean(data5Y))^2))^{1/2}
round(rmse5, 2)
round(colMeans(colMeans(TASMDsummary5)), 2)
GMIM5.result <- colMeans(GMIM5)
round(GMIM5.result, 2)
