# Load necessary libraries
library(readr)         # for reading data files
library(lubridate)     # for handling dates
library(dplyr)         # for data manINDPROulation
library(ggplot2)       # for plotting
library(lmtest)        # for tests (Breusch-Pagan, Durbin-Watson, BG test)
library(sandwich)      # for robust standard errors
library(nlme)          # for gls (ML estimation with variance functions)

#-----------------------------
# Read and prepare the data
#-----------------------------
# Read the CSV file (adjust the path and date format if needed)
data <- read.csv("~/uni/year2/Econometrics 2/USindustrialProduction.csv", stringsAsFactors = FALSE)
data$Date <- as.Date(data$date, format = "%m/%d/%Y")  # adjust the date format as needed

# Create a time trend (assuming the data are ordered by date)
data <- data %>% arrange(Date)
data$Trend <- 1:nrow(data)

# Create a month variable (for seasonal dummies)
data$Month <- factor(month(data$Date))

# Create the log of industrial production (for the log-linear model)
data$logINDPRO <- log(data$INDPRO)

#-----------------------------
# (a) Linear Trend Model
#-----------------------------
# Model: INDPRO = α + β*Trend + ε
model_a <- lm(INDPRO ~ Trend, data = data)
summary(model_a)

# Extract β and its t-statistic for Trend:
beta_a <- coef(summary(model_a))["Trend", "Estimate"]
t_a <- coef(summary(model_a))["Trend", "t value"]
cat("Part (a): β =", beta_a, "with t-statistic =", t_a, "\n")
# Interpretation: β represents the average change in industrial production per time period.

#-----------------------------
# (b) Seasonal dummies without trend
#-----------------------------
# Model: INDPRO = α + seasonal dummies + ε
model_b <- lm(INDPRO ~ Month, data = data)
summary(model_b)

# Report R-squared:
R2_b <- summary(model_b)$r.squared
cat("Part (b): R-squared =", R2_b, "\n")
# Discussion: If the seasonal dummies are highly significant and R² is high,
# the data may not have been seasonally adjusted. (You would elaborate in your answer.)

#-----------------------------
# (c) Linear Trend on Log(INDPRO)
#-----------------------------
# Model: log(INDPRO) = α + β*Trend + ε
model_c <- lm(logINDPRO ~ Trend, data = data)
summary(model_c)

# Extract β and its t-statistic:
beta_c <- coef(summary(model_c))["Trend", "Estimate"]
t_c <- coef(summary(model_c))["Trend", "t value"]
cat("Part (c): β =", beta_c, "with t-statistic =", t_c, "\n")
# Interpretation: Here, β approximates the continuous (percentage) growth rate per period.

#-----------------------------
# (d) Model Preference
#-----------------------------
# Compare the fits (for example, based on R² and residual plots) of the linear (a) versus log-linear (c) model.
# Write up a discussion in your report – e.g., if the log-linear model better captures proportional growth,
# you may prefer it over the linear model.

# For subsequent analysis we continue with the log-linear model.

#-----------------------------
# (e) Regression with Trend and Recession Indicator
#-----------------------------
# Model: log(INDPRO) = α + β1*Trend + β2*RECESSION + ε
model_e <- lm(logINDPRO ~ Trend + RECESSION, data = data)
summary(model_e)
# Report the coefficients and t-statistics for Trend and RECESSION.

#-----------------------------
# (f) Heteroskedasticity Tests
#-----------------------------
# Breusch-Pagan Test:
bp_test <- bptest(model_e)
cat("Part (f) Breusch-Pagan Test:\n")
print(bp_test)

# White Test: Using a variant of BP that includes squared fitted values:
white_test <- bptest(model_e, varformula = ~ fitted(model_e) + I(fitted(model_e)^2))
cat("Part (f) White Test:\n")
print(white_test)

#-----------------------------
# (g) ML Estimation with Variance Model
#-----------------------------
# We assume that the variance is modeled as σi^2 = exp(γ1 + γ2*RECESSION).
# One way is to use the gls() function from nlme with a variance function.
model_g <- gls(logINDPRO ~ Trend, data = data, weights = varExp(form = ~ RECESSION))
summary(model_g)
# In the summary, the parameter associated with the variance function (often labeled as "delta" or similar)
# can be interpreted. γ2 (or its equivalent) and its z-statistic can be extracted from the output.
# (Interpret the significance of this parameter in your report.)

#-----------------------------
# (h) Testing for Serial Correlation
#-----------------------------
# Plot the autocorrelation function (ACF) of the residuals from model (e):
acf(resid(model_e), main = "ACF of Residuals from Model (e)")

# Durbin-Watson Test:
dw <- dwtest(model_e)
cat("Part (h) Durbin-Watson Test:\n")
print(dw)

# Alternatively, use the Breusch-Godfrey test:
bg <- bgtest(model_e, order = 4)
cat("Part (h) Breusch-Godfrey Test:\n")
print(bg)
# In your answer, describe that these tests check for the presence of serial correlation in the residuals.

#-----------------------------
# (i) Including Lagged Dependent Variable
#-----------------------------
# Create a lagged variable for log(INDPRO). Note: The first observation will be NA.
data$lag_logINDPRO <- lag(data$logINDPRO, 1)

# Run the regression including the lagged log(INDPRO):
model_i <- lm(logINDPRO ~ Trend + RECESSION + lag_logINDPRO, data = data)
summary(model_i)

# Extract the coefficient and t-statistic for the lagged dependent variable:
coef_lag <- coef(summary(model_i))["lag_logINDPRO", c("Estimate", "t value")]
cat("Part (i): Lagged log(INDPRO) coefficient =", coef_lag["Estimate"],
    "with t-statistic =", coef_lag["t value"], "\n")
# Explanation: The lagged trend is not added because Trend is a deterministic variable (a time index)
# and including its lag would lead to perfect collinearity or redundant information.

#-----------------------------
# (j) Cochrane-Orcutt Procedure
#-----------------------------
# We start with the model from part (e) (logINDPRO ~ Trend + RECESSION) and then apply Cochrane-Orcutt.
model_e <- lm(logINDPRO ~ Trend + RECESSION, data = data)
summary(model_e)

# Extract residuals
res <- residuals(model_e)

# Estimate the autocorrelation coefficient (rho)
rho_est <- cor(res[-1], res[-length(res)])
cat("Estimated rho =", rho_est, "\n")

# Create a dataset for the transformed model (dropping the first observation)
data_co <- data[-1, ]
# Note: use the corresponding lagged values from the original data
data_co$logINDPRO_lag <- data$logINDPRO[-nrow(data)]
data_co$Trend_lag <- data$Trend[-nrow(data)]
data_co$RECESSION_lag <- data$RECESSION[-nrow(data)]

# Transform variables: subtract rho * lagged value from the current value
data_co$y_star <- data_co$logINDPRO - rho_est * data$logINDPRO[-nrow(data)]
data_co$Trend_star <- data_co$Trend - rho_est * data$Trend[-nrow(data)]
data_co$RECESSION_star <- data_co$RECESSION - rho_est * data$RECESSION[-nrow(data)]

# Re-estimate the model using the transformed data
model_co <- lm(y_star ~ Trend_star + RECESSION_star, data = data_co)
summary(model_co)
# The output will report the estimated autocorrelation parameter (ρ) as well as the adjusted estimates for β.
# In your report, extract γ (i.e. ρ) and the new β estimates and their t-statistics.
