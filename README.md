# US Industrial Production Analysis

A comprehensive analysis of US industrial production trends, seasonality, and recession effects using econometric techniques in R.

## Overview

This project analyzes monthly US industrial production data (from January 1959 to November 2024) to uncover long-term trends, seasonal patterns, and the impact of recessions on production. The analysis uses various econometric models to provide insights that are useful for economists and policymakers.

## Data

The dataset `USindustrialProduction.csv` includes:
- **Date:** The date of the observation (format: MM-DD-YYYY).
- **IP:** The industrial production index.
- **RECESSION:** A binary indicator (0 or 1) showing whether the period is a recession.

Data cleaning steps include converting the date to a proper date format, sorting the data chronologically, creating a time trend variable, and extracting the month.

## Methodology

1. **Trend Analysis:**  
   We first estimate a linear trend model to capture overall changes in production over time.

2. **Seasonality:**  
   We include monthly dummy variables to detect seasonal patterns in the production data.

3. **Log-Linear Modeling:**  
   Transforming the production index with a natural logarithm allows us to interpret changes as percentage growth rates.

4. **Recession Impact:**  
   We add a recession indicator to see how economic downturns affect production levels.

5. **Diagnostic Tests:**  
   We perform tests for heteroskedasticity and serial correlation (e.g., Breusch-Pagan, White, Durbin-Watson) to ensure our model assumptions hold.

6. **Autocorrelation Correction:**  
   We implement a manual Cochrane–Orcutt procedure to adjust for autocorrelation in the error terms.

