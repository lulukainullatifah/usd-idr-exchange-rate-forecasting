# ====================================================================
# USD/IDR EXCHANGE RATE VOLATILITY ANALYSIS
# ====================================================================
# Author  : Luluk Ainul Latifah
# Date    : 2025
# Purpose : Forecast USD/IDR volatility using ARMA-TARCH & ARMA-EGARCH
# Data    : Daily USD/IDR rates, Bank Indonesia (Jan 2021 - Dec 2024)
# ====================================================================

# --------------------------------------------------------------------
# PART A: SETUP — Load Required Packages
# --------------------------------------------------------------------
library(nloptr)   # For Genetic Algorithm solver (nloptr)
library(tseries)  # For ADF test and Jarque-Bera test
library(e1071)    # For skewness() and kurtosis()
library(FinTS)    # For ArchTest() — ARCH effect testing
library(rugarch)  # For TARCH and EGARCH modeling
library(ggplot2)  # For custom visualizations
library(lmtest)   # For coeftest() — clean coefficient tables
library(xts)      # For time-indexed log return series

# NOTE: Data must be loaded as dataframe 'data2' with columns:
#   - Tanggal : Date column (Date format)
#   - Kurs    : USD/IDR exchange rate (numeric)

# --------------------------------------------------------------------
# PART B: DATA DESCRIPTION AND TRANSFORMATION
# --------------------------------------------------------------------

# --- B1. Descriptive Statistics ---
cat("=== PART B: DATA DESCRIPTION ===\n")
print(summary(data2$Kurs))
cat("Standard Deviation :", sd(data2$Kurs, na.rm = TRUE), "\n")
cat("Skewness           :", skewness(data2$Kurs, na.rm = TRUE), "\n")
cat("Kurtosis           :", kurtosis(data2$Kurs, na.rm = TRUE), "\n")
cat("\nJarque-Bera Normality Test:\n")
print(jarque.bera.test(data2$Kurs))

# --- B2. Raw Data Visualization ---
plot(data2$Tanggal, data2$Kurs,
     type = 'l',
     main = "Daily USD/IDR Exchange Rate (2021-2024)",
     xlab = "Year",
     ylab = "Exchange Rate (IDR)")

hist(data2$Kurs,
     main  = "Histogram of USD/IDR Exchange Rate",
     xlab  = "Exchange Rate (IDR)",
     ylab  = "Frequency",
     col   = "lightblue",
     border = "black")

# --- B3. Stationarity Test on Raw Data ---
cat("\n=== ADF Test on Raw Data ===\n")
print(adf.test(data2$Kurs))
# Expected: p-value > 0.05 → non-stationary → transformation required

# --- B4. Log Return Transformation ---
# Formula: r_t = ln(P_t / P_{t-1})
data2$log_return  <- c(NA, diff(log(data2$Kurs)))
log_return_clean  <- na.omit(data2$log_return)

# Create xts object with date index for rugarch compatibility
log_return_xts <- xts(log_return_clean,
                      order.by = data2$Tanggal[!is.na(data2$log_return)])

# --- B5. Stationarity Test on Log Returns ---
cat("\n=== ADF Test on Log Returns ===\n")
print(adf.test(log_return_clean))
# Expected: p-value < 0.05 → stationary ✓

plot(log_return_clean,
     type = "l",
     main = "Log Return of USD/IDR Exchange Rate",
     ylab = "Log Return",
     xlab = "Observation Period")

# --------------------------------------------------------------------
# PART C: MEAN MODEL — ARMA IDENTIFICATION AND SELECTION
# --------------------------------------------------------------------
cat("\n=== PART C: ARMA MODEL IDENTIFICATION ===\n")

# --- C1. ACF and PACF Plots ---
par(mfrow = c(1, 2))
acf(log_return_clean,  lag.max = 36, main = "ACF of Log Returns")
pacf(log_return_clean, lag.max = 36, main = "PACF of Log Returns")
par(mfrow = c(1, 1))

# --- C2. Candidate Model Comparison ---
model_candidates <- list(
  "AR(1)"     = c(1, 0, 0),
  "MA(1)"     = c(0, 0, 1),
  "ARMA(1,1)" = c(1, 0, 1),
  "ARMA(1,2)" = c(1, 0, 2),
  "ARMA(2,1)" = c(2, 0, 1),
  "ARMA(2,2)" = c(2, 0, 2)
)

for (model_name in names(model_candidates)) {
  cat("\n---", model_name, "---\n")
  fit <- arima(log_return_clean, order = model_candidates[[model_name]])
  print(coeftest(fit))
  cat("AIC:", AIC(fit), "\n")
  cat("BIC:", BIC(fit), "\n")
  cat("Ljung-Box Test:\n")
  print(Box.test(fit$residuals, type = "Ljung-Box"))
  cat("ARCH-LM Test:\n")
  print(ArchTest(fit$residuals))
}

# --- C3. Best ARMA Model: ARMA(2,2) ---
# Selected based on: lowest AIC (-8397.90), all parameters significant,
# residuals pass white noise and ARCH-LM tests
model_arma_best  <- arima(log_return_clean, order = c(2, 0, 2))
residuals_arma   <- model_arma_best$residuals
cat("\n✓ ARMA(2,2) selected as best mean model.\n")

# ACF & PACF of Squared Residuals (check for ARCH effects)
par(mfrow = c(1, 2))
acf(residuals_arma^2,  lag.max = 36, main = "ACF of Squared Residuals")
pacf(residuals_arma^2, lag.max = 36, main = "PACF of Squared Residuals")
par(mfrow = c(1, 1))

# --------------------------------------------------------------------
# PART D: ASYMMETRY TEST — Sign Bias Test
# --------------------------------------------------------------------
cat("\n=== PART D: SIGN BIAS TEST ===\n")

spec_test <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE)
)
fit_test <- ugarchfit(spec = spec_test, data = residuals_arma, solver = 'hybrid')

if (convergence(fit_test) == 0) {
  show(fit_test)
  cat("\nINTERPRETATION:\n")
  cat("  If 'Negative Sign Bias' or 'Joint Effect' p-value < 0.05\n")
  cat("  → Asymmetric effects confirmed → use TARCH or EGARCH\n")
} else {
  cat("WARNING: Sign Bias test model failed to converge.\n")
}
# Result: Negative Sign Bias p=0.0017, Joint Effect p=0.0079
# → Asymmetric effects CONFIRMED → TARCH & EGARCH appropriate ✓

# --------------------------------------------------------------------
# PART E: VOLATILITY MODEL SPECIFICATIONS
# --------------------------------------------------------------------

# --- E1. Model Specifications ---
spec_tarch <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(2, 1)),
  mean.model         = list(armaOrder = c(2, 2), include.mean = TRUE),
  distribution.model = "std"   # Student's t distribution for fat tails
)

spec_egarch <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(2, 1)),
  mean.model         = list(armaOrder = c(2, 2), include.mean = TRUE),
  distribution.model = "std"
)

# --- E2. Optimization Algorithm Comparison ---
cat("\n=== PART E: OPTIMIZATION ALGORITHM COMPARISON ===\n")

model_specs <- list(
  "ARMA(2,2)-TARCH(2,1)"  = spec_tarch,
  "ARMA(2,2)-EGARCH(2,1)" = spec_egarch
)

solvers_to_test <- list(
  "BFGS (hybrid)"             = list(name = "hybrid", control = list()),
  "Newton-Raphson (solnp)"    = list(name = "solnp",  control = list()),
  "Genetic Algorithm (nloptr)"= list(name = "nloptr",
                                     control = list(solver   = "NLOPT_GN_ISRES",
                                                    xtol_rel = 1e-8,
                                                    maxeval  = 10000))
)

hasil_optimasi <- data.frame(
  Model         = character(),
  Solver        = character(),
  LogLikelihood = numeric(),
  Convergence   = character(),
  Time_sec      = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in names(model_specs)) {
  for (solver_name in names(solvers_to_test)) {
    cat(paste("\n--> Estimating:", model_name, "| Solver:", solver_name, "\n"))

    current_spec    <- model_specs[[model_name]]
    current_solver  <- solvers_to_test[[solver_name]]$name
    current_control <- solvers_to_test[[solver_name]]$control

    start_time   <- Sys.time()
    fit_attempt  <- try({
      ugarchfit(spec           = current_spec,
                data           = log_return_clean,
                solver         = current_solver,
                solver.control = current_control)
    }, silent = TRUE)
    end_time       <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (!inherits(fit_attempt, "try-error")) {
      loglik      <- fit_attempt@fit$LLH
      conv_code   <- convergence(fit_attempt)
      conv_status <- ifelse(conv_code == 0, "Converged", "Failed")
    } else {
      loglik      <- NA
      conv_status <- "Error"
    }

    hasil_optimasi <- rbind(hasil_optimasi, data.frame(
      Model         = model_name,
      Solver        = solver_name,
      LogLikelihood = loglik,
      Convergence   = conv_status,
      Time_sec      = round(execution_time, 2)
    ))
    cat(paste("   Status:", conv_status,
              "| Log-Lik:", round(loglik, 3),
              "| Time:", round(execution_time, 2), "sec\n"))
  }
}

cat("\n=== OPTIMIZATION COMPARISON RESULTS ===\n")
print(hasil_optimasi[order(hasil_optimasi$Model, -hasil_optimasi$LogLikelihood), ])
# Result: BFGS and Newton-Raphson converge successfully
# Genetic Algorithm fails due to log-likelihood surface complexity
# → BFGS (hybrid) selected as estimation method

# --------------------------------------------------------------------
# PART F: PARAMETER ESTIMATION AND DIAGNOSTIC VALIDATION
# --------------------------------------------------------------------
solver_terpilih <- 'hybrid'

# --- F1. TARCH(2,1) Estimation ---
cat("\n=== PART F1: ARMA(2,2)-TARCH(2,1) ===\n")
fit_tarch <- ugarchfit(spec   = spec_tarch,
                       data   = log_return_xts,
                       solver = solver_terpilih)
show(fit_tarch)

if (convergence(fit_tarch) == 0) {
  par(mfrow = c(2, 2))
  plot(fit_tarch, which = 1)   # Returns with Conditional SD
  plot(fit_tarch, which = 8)   # Empirical Density
  plot(fit_tarch, which = 9)   # QQ-Plot
  plot(fit_tarch, which = 7)   # Cross-Correlation of Squared Residuals
  par(mfrow = c(1, 1))
}
# Note: TARCH passes all diagnostic tests BUT parameters largely
# insignificant + unstable (Nyblom Joint Stat = 315.89 >> critical 2.96)

# --- F2. EGARCH(2,1) Estimation ---
cat("\n=== PART F2: ARMA(2,2)-EGARCH(2,1) ===\n")
fit_egarch <- ugarchfit(spec   = spec_egarch,
                        data   = log_return_xts,
                        solver = solver_terpilih)
show(fit_egarch)

if (convergence(fit_egarch) == 0) {
  par(mfrow = c(2, 2))
  plot(fit_egarch, which = 1)
  plot(fit_egarch, which = 8)
  plot(fit_egarch, which = 9)
  plot(fit_egarch, which = 7)
  par(mfrow = c(1, 1))
}
# EGARCH results: Parameters stable (Nyblom = 2.55 < 2.96)
# alpha1a (leverage) significant → confirms asymmetric effect
# beta1 = 0.970 → high volatility persistence ✓

# --------------------------------------------------------------------
# PART G: FORECASTING AND EVALUATION
# --------------------------------------------------------------------

# --- G1. Out-of-Sample MSE Evaluation ---
cat("\n=== PART G1: MSE EVALUATION (Out-of-Sample) ===\n")

n_eval  <- 10
n_total <- length(log_return_clean)
train_data <- log_return_clean[1:(n_total - n_eval)]
test_data  <- log_return_clean[(n_total - n_eval + 1):n_total]

fit_tarch_train  <- ugarchfit(spec = spec_tarch,  data = train_data, solver = solver_terpilih)
fit_egarch_train <- ugarchfit(spec = spec_egarch, data = train_data, solver = solver_terpilih)

forecast_tarch_eval  <- ugarchforecast(fit_tarch_train,  n.ahead = n_eval)
forecast_egarch_eval <- ugarchforecast(fit_egarch_train, n.ahead = n_eval)

mse_tarch  <- mean((test_data - fitted(forecast_tarch_eval))^2)
mse_egarch <- mean((test_data - fitted(forecast_egarch_eval))^2)

hasil_mse <- data.frame(
  Model              = c("ARMA(2,2)-TARCH(2,1)", "ARMA(2,2)-EGARCH(2,1)"),
  Mean_Squared_Error = c(mse_tarch, mse_egarch)
)
cat("\nMSE Comparison:\n")
print(hasil_mse)
# TARCH MSE  = 2.343820e-05 (slightly lower)
# EGARCH MSE = 2.366173e-05
# → EGARCH preferred despite marginally higher MSE due to robustness

# --- G2. 5-Day Ahead Forecast ---
cat("\n=== PART G2: 5-DAY AHEAD FORECAST (Jan 2-8, 2025) ===\n")

n_horizon <- 5
forecast_tarch_final  <- ugarchforecast(fit_tarch,  n.ahead = n_horizon)
forecast_egarch_final <- ugarchforecast(fit_egarch, n.ahead = n_horizon)

# Extract forecast components
return_tarch  <- as.numeric(fitted(forecast_tarch_final))
vol_tarch     <- as.numeric(sigma(forecast_tarch_final))
return_egarch <- as.numeric(fitted(forecast_egarch_final))
vol_egarch    <- as.numeric(sigma(forecast_egarch_final))

# 95% Confidence Interval (using Normal quantiles)
alpha <- 0.05
z_low <- qnorm(alpha / 2)        # ≈ -1.96
z_up  <- qnorm(1 - alpha / 2)    # ≈ +1.96

lower_tarch  <- return_tarch  + z_low * vol_tarch
upper_tarch  <- return_tarch  + z_up  * vol_tarch
lower_egarch <- return_egarch + z_low * vol_egarch
upper_egarch <- return_egarch + z_up  * vol_egarch

# Forecast dates (business days, excluding Jan 1 holiday)
tanggal_final <- as.Date(c(
  "2025-01-02", "2025-01-03", "2025-01-06",
  "2025-01-07", "2025-01-08"
))

# Compile forecast table
tabel_peramalan_final <- data.frame(
  Period                 = 1:n_horizon,
  Date                   = tanggal_final,
  TARCH_Return_Forecast  = return_tarch,
  TARCH_Lower_95CI       = lower_tarch,
  TARCH_Upper_95CI       = upper_tarch,
  TARCH_Volatility       = vol_tarch,
  EGARCH_Return_Forecast = return_egarch,
  EGARCH_Lower_95CI      = lower_egarch,
  EGARCH_Upper_95CI      = upper_egarch,
  EGARCH_Volatility      = vol_egarch
)

cat("\n5-Day Ahead Forecast Results:\n")
print(tabel_peramalan_final, row.names = FALSE)

# Back-transform log returns to IDR exchange rate levels
# Formula: P_t = P_{t-1} * exp(r_t)
last_rate <- tail(data2$Kurs, 1)

kurs_tarch  <- numeric(n_horizon)
kurs_egarch <- numeric(n_horizon)

for (i in 1:n_horizon) {
  prev_rate     <- if (i == 1) last_rate else kurs_tarch[i - 1]
  kurs_tarch[i] <- prev_rate * exp(return_tarch[i])
}
for (i in 1:n_horizon) {
  prev_rate      <- if (i == 1) last_rate else kurs_egarch[i - 1]
  kurs_egarch[i] <- prev_rate * exp(return_egarch[i])
}

tabel_kurs_forecast <- data.frame(
  Date             = tanggal_final,
  TARCH_Rate_IDR   = round(kurs_tarch,  2),
  EGARCH_Rate_IDR  = round(kurs_egarch, 2)
)

cat("\nForecasted USD/IDR Exchange Rate (IDR):\n")
print(tabel_kurs_forecast, row.names = FALSE)
cat("\n=== ANALYSIS COMPLETE ===\n")
