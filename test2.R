library(readstata13)

# ----------------------------------------------------------
# Comprehensive Test File for FDID Package
# ----------------------------------------------------------

# Set a seed for reproducibility
set.seed(12345)

# Load the FDID package and required libraries
library(fdid)
library(haven)
library(dplyr)

d <- read.dta13("Table1.dta")
wubaos <- names(tapply(d$stronghold_reform, d$prefecture_id, mean)
                [tapply(d$stronghold_reform, d$prefecture_id, mean)>0])
d$wubao <- ifelse(d$prefecture_id %in% wubaos, 1, 0)

d$uniqueid <- paste(as.character(d$prov_id), as.character(d$prefecture_id), sep = "-")
d$prefecture_id = factor(d$prefecture_id)
# Define common parameters for all tests
start <- 2
end <- 12
period <- start:end
tr_period <- 8:12
treat <- "wubao"  # Treatment variable
Y <- "arist_recruit"     # Outcome variable
time <- "emperor_id"
unit <- "uniqueid"
cluster <- "prefecture_id"

# ------------------
# 1. Test: fdid_prepare()
# ------------------

# Basic test: Prepare data and verify structure
prepared_data <- fdid_prepare(
  data = d,
  Y_label = Y,
  X_labels = NULL,
  G_label = treat,
  unit_label = unit,
  time_label = time,
  cluster_label = cluster
)
fdid_results_ols1 <- fdid(
  s = prepared_data,
  tr_period = tr_period,
  ref_period = 2:7,
  method = "ols1",
  vartype = "robust"
)

plot(fdid_results_ols1)
plot(fdid_results_ols1, type="dynamic")
