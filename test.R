# ----------------------------------------------------------
# Comprehensive Test File for FDID Package
# ----------------------------------------------------------

# Set a seed for reproducibility

# Load the FDID package and required libraries
library(fdid)
library(haven)
library(dplyr)

# Load the dataset
d <- as.data.frame(read.csv("mortality.csv"))

# Data Preparation: Add new variables and identifiers
d$uniqueid <- paste(as.character(d$provid), as.character(d$countyid), sep = "-")
d$provid <- as.factor(d$provid)
d$highzupu <- ifelse(d$pczupu >= median(d$pczupu, na.rm = TRUE), 1, 0)

# Define common parameters for all tests
start <- 1954
end <- 1966
period <- start:end
tr_period <- 1958:1961
treat <- "highzupu"  # Treatment variable
Y <- "mortality"     # Outcome variable
time <- "year"
unit <- "uniqueid"
cluster <- "provid"
covar <- c("avggrain", "nograin", "urban", "dis_bj", "dis_pc",
           "rice", "minority", "edu", "lnpop")

# ------------------
# 1. Test: fdid_prepare()
# ------------------

# Basic test: Prepare data and verify structure
prepared_data <- fdid_prepare(
  data = d,
  Y_label = Y,
  X_labels = covar,
  G_label = treat,
  unit_label = unit,
  time_label = time,
  # cluster_label = cluster
)
cat("Prepared Data Summary:\n")
print(str(prepared_data))

# Edge case: Missing values in key columns
d_missing <- d
d_missing$avggrain[1:10] <- NA
prepared_data_missing <- fdid_prepare(
  data = d_missing,
  Y_label = Y,
  X_labels = covar,
  G_label = treat,
  unit_label = unit,
  time_label = time
)
cat("Prepared Data with Missing Values:\n")
print(head(prepared_data_missing))

# ------------------
# 2. Test: fdid()
# ------------------

# Test 2.1: Basic FDID estimation with OLS1
fdid_results_ols1 <- fdid(
  s = prepared_data,
  tr_period = tr_period,
  ref_period = 1954:1957,
  method = "ols1",
  vartype = "robust"
)
cat("FDID Results Summary (OLS1):\n")
summary(fdid_results_ols1)
cat("FDID Object Overview (OLS1):\n")
print(fdid_results_ols1)

# Test 2.2: FDID with Entropy Balancing
fdid_results_ebal <- fdid(
  s = prepared_data,
  tr_period = tr_period,
  ref_period = 1954:1957,
  method = "ebal",
  vartype = "robust"
)
cat("FDID Results Summary (Entropy Balancing):\n")
summary(fdid_results_ebal)

# Test 2.3: FDID with AIPW
fdid_results_aipw <- fdid(
  s = prepared_data,
  tr_period = tr_period,
  ref_period = 1954:1957,
  method = "aipw",
  vartype = "robust"
)
cat("FDID Results Summary (AIPW):\n")
summary(fdid_results_aipw)

# Test 2.4: Missing values handling in FDID
fdid_results_missing <- fdid(
  s = prepared_data_missing,
  tr_period = tr_period,
  ref_period = 1954:1957,
  method = "ols1",
  vartype = "robust",
  na.rm = TRUE
)
cat("FDID Results Summary with Missing Values (OLS1):\n")
summary(fdid_results_missing)

# Test 2.5: Parallel computation with bootstrap variance
fdid_results_parallel <- fdid(
  s = prepared_data,
  tr_period = tr_period,
  ref_period = 1954:1957,
  method = "ols1",
  vartype = "bootstrap",
  parallel = TRUE,
  cores = 8,
  nsims = 1000
)
cat("FDID Results Summary (Parallel Bootstrap):\n")
summary(fdid_results_parallel)

# ------------------
# 3. Test: plot.fdid()
# ------------------

# Test 3.1: Raw means plot
plot.fdid(fdid_results_ols1, type = "raw")

# Test 3.2: Dynamic effects plot
plot.fdid(fdid_results_ols1, type = "dynamic")

# Test 3.3: Propensity score overlap
plot.fdid(fdid_results_aipw, type = "overlap")

# Test 3.4: Comparison of estimates
fdid_results_ols2 <- fdid(
  s = prepared_data,
  tr_period = tr_period,
  pre_period = 1954:1957,
  method = "ols2",
  vartype = "robust"
)
plot.fdid(list(fdid_results_ols1, fdid_results_ols2), type = "est")

# ------------------
# 4. Test: Edge Cases
# ------------------

# Edge Case 1: No covariates provided
prepared_data_nocovar <- fdid_prepare(
  data = d,
  Y_label = Y,
  X_labels = character(0),  # No covariates
  G_label = treat,
  unit_label = unit,
  time_label = time
)
fdid_results_nocovar <- fdid(
  s = prepared_data_nocovar,
  tr_period = tr_period,
  ref_period = 1954:1957,
  method = "ols1",
  vartype = "robust"
)
cat("FDID Results Summary (No Covariates):\n")
summary(fdid_results_nocovar)

# Edge Case 2: Single treatment period
fdid_results_single <- fdid(
  s = prepared_data,
  tr_period = c(1958),
  ref_period = c(1957),
  method = "ols1",
  vartype = "robust"
)
cat("FDID Results Summary (Single Treatment Period):\n")
summary(fdid_results_single)

