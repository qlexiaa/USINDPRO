install.packages("sandwich")
# Load necessary libraries
library(tidyverse)  # For data manipulation
library(lmtest)     # For Breusch-Pagan and White tests
library(sandwich)   # For robust standard errors (if needed)
library(stats)      # For basic statistical functions

# Read the data
data <- read.csv("~/uni/year2/Econometrics 2/USindustrialProduction.csv", stringsAsFactors = FALSE)

# Convert date to "day month year" format (assuming input is still MM/DD/YYYY)
# First, parse the original MM/DD/YYYY format
data$date <- as.Date(data$date, format = "%m/%d/%Y")

# Then, reformat to "DD Month YYYY" (e.g., "01 January 1959")
data$date_formatted <- format(data$date, "%d %B %Y")

# Create a time trend variable (i = 1, 2, ..., n)
data$trend <- 1:nrow(data)

# (a) Linear trend model: y_i = alpha + beta * i + epsilon_i
model_a <- lm(INDPRO ~ trend, data = data)
summary_a <- summary(model_a)
beta_a <- coef(model_a)["trend"]
t_stat_a <- summary_a$coefficients["trend", "t value"]
cat("Question (a):\n")
cat("Beta (trend coefficient):", beta_a, "\n")
cat("t-statistic:", t_stat_a, "\n")
# Interpretation: Beta represents the average monthly change in industrial production.

# (b) Seasonal dummies model (no trend)
data$month <- format(data$date, "%m")  # Extract month (still numeric for dummies)
model_b <- lm(INDPRO ~ factor(month), data = data)
r_squared_b <- summary(model_b)$r.squared
cat("\nQuestion (b):\n")
cat("R-squared:", r_squared_b, "\n")
# Check if data is seasonally adjusted
seasonal_variation <- anova(model_b)$"Pr(>F)"[1] < 0.05
cat("Was data seasonally adjusted? ", ifelse(!seasonal_variation, "Likely yes", "Likely no"), 
    " (based on weak seasonal effects if p-value > 0.05).\n")

# (c) Linear trend model with log of industrial production
data$log_INDPRO <- log(data$INDPRO)
model_c <- lm(log_INDPRO ~ trend, data = data)
summary_c <- summary(model_c)
beta_c <- coef(model_c)["trend"]
t_stat_c <- summary_c$coefficients["trend", "t value"]
cat("\nQuestion (c):\n")
cat("Beta (trend coefficient):", beta_c, "\n")
cat("t-statistic:", t_stat_c, "\n")
# Interpretation: Beta is the approximate monthly percentage change in industrial production.

# (d) Model preference (linear vs log-linear)
cat("\nQuestion (d):\n")
cat("Preference: Compare R^2 and residuals. Log-linear often better for trending data.\n")
cat("Linear R^2:", summary_a$r.squared, "vs Log-linear R^2:", summary_c$r.squared, "\n")
cat("Use diagnostic plots (e.g., qqnorm) to confirm. Log-linear typically preferred for economic time series.\n")

# Generate QQ plots for both models
par(mfrow = c(1, 2))  # Set up a 1x2 plotting grid
qqnorm(residuals(model_a), main = "QQ Plot: Linear Model (INDPRO ~ trend)")
qqline(residuals(model_a), col = "red")
qqnorm(residuals(model_c), main = "QQ Plot: Log-Linear Model (log_INDPRO ~ trend)")
qqline(residuals(model_c), col = "red")
par(mfrow = c(1, 1))  # Reset plotting grid
cat("QQ Plot Interpretation: Check if points follow the red line. Closer alignment suggests better normality of residuals.\n")


# First differences of log_INDPRO
data$diff_log_INDPRO <- c(NA, diff(data$log_INDPRO))

# (e) Log-linear model with trend and recession indicator
model_e <- lm(log_INDPRO ~ trend + RECESSION, data = data)
summary_e <- summary(model_e)
beta_trend_e <- coef(model_e)["trend"]
t_stat_trend_e <- summary_e$coefficients["trend", "t value"]
beta_recession_e <- coef(model_e)["RECESSION"]
t_stat_recession_e <- summary_e$coefficients["RECESSION", "t value"]
cat("\nQuestion (e):\n")
cat("Trend coefficient:", beta_trend_e, "t-statistic:", t_stat_trend_e, "\n")
cat("Recession coefficient:", beta_recession_e, "t-statistic:", t_stat_recession_e, "\n")

# (f) Test for heteroskedasticity
bp_test <- bptest(model_e)
white_test <- bptest(model_e, ~ fitted(model_e) + I(fitted(model_e)^2))
cat("\nQuestion (f):\n")
cat("Breusch-Pagan test statistic:", bp_test$statistic, "p-value:", bp_test$p.value, "\n")
cat("White test statistic:", white_test$statistic, "p-value:", white_test$p.value, "\n")

# (g) ML estimation with heteroskedasticity: sigma^2 = exp(gamma1 + gamma2 * recession)
loglik <- function(par, y, X, recession) {
  alpha <- par[1]
  beta <- par[2]
  gamma1 <- par[3]
  gamma2 <- par[4]
  mu <- X %*% c(alpha, beta)
  sigma2 <- exp(gamma1 + gamma2 * recession)
  -sum(dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE))
}
X <- model.matrix(~ trend, data = data)
init_par <- c(coef(model_e)[1], coef(model_e)[2], 0, 0)  # Initial values
ml_fit <- optim(init_par, loglik, y = data$log_INDPRO, X = X, recession = data$RECESSION, 
                method = "BFGS", hessian = TRUE)
gamma2 <- ml_fit$par[4]
se_gamma2 <- sqrt(diag(solve(ml_fit$hessian)))[4]
z_stat_gamma2 <- gamma2 / se_gamma2
cat("\nQuestion (g):\n")
cat("Gamma2:", gamma2, "z-statistic:", z_stat_gamma2, "\n")
cat("Significant if |z| > 1.96 (5% level).\n")

# (h) Test for serial correlation
acf_resid <- acf(residuals(model_e), main = "ACF of Residuals")
lag_length <- max(which(abs(acf_resid$acf[-1]) > 2/sqrt(nrow(data))))  # Significant lags
dw_test <- dwtest(model_e)
cat("\nQuestion (h):\n")
cat("Lag length from ACF:", lag_length, "\n")
cat("Durbin-Watson test statistic:", dw_test$statistic, "p-value:", dw_test$p.value, "\n")
cat("Motivation: DW test checks AR(1) serial correlation, suitable for time series.\n")

# (i) Add lagged dependent variable
data$lag_log_INDPRO <- lag(data$log_INDPRO)
model_i <- lm(log_INDPRO ~ trend + RECESSION + lag_log_INDPRO, data = data)
summary_i <- summary(model_i)
coef_lag <- coef(model_i)["lag_log_INDPRO"]
t_stat_lag <- summary_i$coefficients["lag_log_INDPRO", "t value"]
cat("\nQuestion (i):\n")
cat("Lagged dependent variable coefficient:", coef_lag, "t-statistic:", t_stat_lag, "\n")
cat("Lagged trend not added: Perfect collinearity with trend (trend_t = trend_{t-1} + 1).\n")

# (j) Cochrane-Orcutt procedure
resid_e <- residuals(model_e)
ar1_model <- lm(resid_e[-1] ~ resid_e[-length(resid_e)])
rho <- coef(ar1_model)[2]
data$y_transformed <- data$log_INDPRO - rho * lag(data$log_INDPRO)
data$trend_transformed <- data$trend - rho * lag(data$trend)
data$recession_transformed <- data$RECESSION - rho * lag(data$RECESSION)
model_j <- lm(y_transformed ~ trend_transformed + recession_transformed, data = data[-1, ])
summary_j <- summary(model_j)
beta_j <- coef(model_j)["trend_transformed"]
t_stat_j <- summary_j$coefficients["trend_transformed", "t value"]
cat("\nQuestion (j):\n")
cat("Rho (AR(1) coefficient):", rho, "\n")
cat("New Beta (trend):", beta_j, "t-statistic:", t_stat_j, "\n")

# Optional: View the first few rows with the new date format
head(data[, c("date_formatted", "INDPRO", "RECESSION")])