# Load data
install.packages("forecast")

data <- read.csv("C:\\Users\\91830\\Downloads\\Telegram Desktop\\monthly_milk_production.csv")

date <- as.vector(data$Date)
qnt <- as.vector(data$Production)
qnt <- qnt[c(1:168)]  # first 168 obs
qnt

# Summary
summary(qnt)

# Time series object
as <- ts(qnt, freq = 12, start = c(1962, 1))
n <- length(qnt)

plot.ts(as, main = "Milk Production Time Series")

###########Test for randomness (Turning Point Test)
u <- 0
for (i in 2:(n - 1)) {
  if (((qnt[i] > qnt[i + 1]) && (qnt[i] > qnt[i - 1])) ||
      ((qnt[i] < qnt[i + 1]) && (qnt[i] < qnt[i - 1]))) {
    u <- u + 1
  }
}
u
e_u <- 2 * (n - 2) / 3
v_u <- (16 * n - 29) / 90
z <- (u - e_u) / sqrt(v_u)
abs(z)

#######333 Test for trend (Relative ordering test) #########
q <- 0
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    if (qnt[i] > qnt[j]) {
      q <- q + 1
    }
  }
}
q
Expecq <- n * (n - 1) / 4
tu <- 1 - 4 * q / (n * (n - 1))
V_tu <- 2 * (2 * n + 5) / (9 * n * (n - 1))
Z <- tu / sqrt(V_tu)
Z  # test statistic

######3333 Seasonality test (Friedman’s Test) #########3
qnt_array <- array(qnt, dim = c(12, 14))
r <- array(0, dim = c(12, 14))
for (j in 1:14) {
  y <- rank(qnt_array[, j])
  r[, j] <- y
}
me <- rowSums(r)
c1 <- 14 * (13) / 2
c2 <- 14 * 13
chi <- sum((me - c1)^2) / c2
chi

#### Differencing for seasonality 
qnt_r <- qnt[13:168] - qnt[1:156]
qnt_array <- array(qnt_r, dim = c(12, 13))
r <- array(0, dim = c(12, 13))
for (j in 1:13) {
  y <- rank(qnt_array[, j])
  r[, j] <- y
}
me2 <- rowSums(r)
c1 <- 13 * (13) / 2
c2 <- 13 * 13
chi <- sum((me2 - c1)^2) / c2
chi

acf(qnt, lag.max = 48)

#### Estimating trend component (Moving Average)
m <- rep(0, 156)
for (i in 7:162) {
  m[i - 6] <- (1 / 12) * (0.5 * qnt[i - 6] + sum(qnt[(i - 5):(i + 5)]) +
                            0.5 * qnt[i + 6])
}
a1 <- ts(m, freq = 12, start = c(1962, 6))
plot.ts(a1, type = 'o', main = "Trend Component", ylab = "estimated trend")

# Removing trend
xr <- rep(0, 156)
for (i in 7:162) {
  xr[i - 6] <- qnt[i] - m[i - 6]
}
a2 <- ts(xr, freq = 12, start = c(1962, 1))
plot.ts(a2, type = 'o', main = "After Removing Trend")

#### Estimate seasonality component 
rmv_sesnlity <- function(x) {
  x <- matrix(x, ncol = 12, byrow = TRUE)
  w <- colMeans(x)
  return(w)
}
w <- rmv_sesnlity(xr)

est_sesnlity <- function(x) {
  s <- x - mean(x)
  return(s)
}
estmtd_sesnlity <- est_sesnlity(w)
plot(estmtd_sesnlity, type = 'o', col = "blue",
     main = "Estimated Seasonality Component")

#### Deseasonalize data
deseasonlzation <- function(x, y) {
  z <- matrix(x, ncol = 12, byrow = TRUE) -
    matrix(rep(y, length(x) / 12), ncol = 12, byrow = TRUE)
  return(as.vector(t(z)))
}
desesnlizd_detrndd_dt <- deseasonlzation(xr, estmtd_sesnlity)

library(forecast)
plot(desesnlizd_detrndd_dt, type = 'o', col = "blue",
     main = "Deseasonalized & Detrended Data")

#### Stationarity test 
if (!require(tseries)) {
  install.packages("tseries")
  library(tseries)
}
adf.test(desesnlizd_detrndd_dt)

#### White noise test
acf(desesnlizd_detrndd_dt, lag.max = 100)

### Decomposition
d <- decompose(as)
plot(d)

y <- na.omit(d$random)
acf(y, lag.max = 110, main = "ACF of Random Component")
pacf(y, lag.max = 100, main = "PACF of Random Component")

### AIC calculation for ARMA models
aic <- matrix(0, 10, 10)
for (p in 0:9) {
  for (q in 0:9) {
    a.p.q <- arima(y, order = c(p, 0, q), method = "ML")
    aic[p + 1, q + 1] <- a.p.q$aic
  }
}
aic

### Final ARMA model 
model <- arima(y, order = c(7, 0, 6), method = "ML")
model$aic
model$coef

### Check stationarity of AR polynomial 
z <- c(1, 0.15434185, -0.25690396, 0.21989901,
       0.23248325, 0.41171631, 0.24343836, 0.10272174)
polyroot(z)

## Forecast
pred <- predict(model, n.ahead = 12, se.fit = TRUE)
t <- na.omit(d$trend)
deep <- tail(t, 12)
mat <- cbind(deep, estmtd_sesnlity, pred$pred)
predicted <- rowSums(mat)

# Compare with actual
actual <- qnt[151:162]
last <- cbind(predicted, actual)
print(last)

# Plot forecast vs actual
x1 <- 1:12
plot(x1, actual, type = "o", col = "red", main = "Forecast vs Actual")
lines(x1, predicted, col = "green")
legend("topleft", legend = c("Actual", "Predicted"),
       col = c("red", "green"), lty = 1)
