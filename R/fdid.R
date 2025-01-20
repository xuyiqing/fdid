#' Fixed Difference-in-Differences Estimation
#'
#' Performs fixed difference-in-differences (FDID) estimation using various methods and variance estimation techniques.
#'
#' @param s A data frame prepared using `fdid_prepare`.
#' @param tr_period A numeric vector specifying the treatment periods.
#' @param ref_period A numeric vector specifying the reference period(s).
#' @param method A string specifying the estimation method. Options: `"ols1"`, `"ols2"`, `"did"`, `"ebal"`, `"aipw"`. Default is `"ols1"`.
#' @param vartype A string specifying the variance estimation type. Options: `"robust"`, `"bootstrap"`, `"jackknife"`. Default is `"robust"`.
#' @param na.rm Logical; whether to remove rows with missing values in relevant columns. Default is `FALSE`.
#' @param nsims Number of simulations for bootstrap variance estimation. Default is `1000`.
#' @param parallel Logical; whether to perform parallel computations. Default is `FALSE`.
#' @param cores Number of cores for parallel computations. Default is `2`.
#'
#' @return A list with the following components:
#' - `static`: Static FDID estimate.
#' - `dynamic`: Dynamic FDID estimates for each year.
#' - `raw_means`: Raw mean outcomes by group and year.
#' - Metadata, including periods, method, variance type, and propensity scores (if applicable).
#'
#' @examples
#' fdid_results <- fdid(
#'   s = prepared_data,
#'   tr_period = c(1958, 1959, 1960, 1961),
#'   ref_period = 1954:1957,
#'   method = "ebal",
#'   vartype = "robust"
#' )
#' print(fdid_results$dynamic)
#'
#' @import dplyr
#' @import estimatr
#' @importFrom grf causal_forest
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
#' @export
fdid <- function(s,
                 tr_period,
                 ref_period,
                 method   = "ols1",
                 vartype  = "robust",
                 na.rm    = FALSE,
                 nsims    = 1000,
                 parallel = FALSE,
                 cores    = 2) {
  # --- Needed packages ---
  if (!requireNamespace("estimatr", quietly=TRUE)) {
    stop("Package 'estimatr' is required but not installed.")
  }
  if (!requireNamespace("car", quietly=TRUE)) {
    stop("Package 'car' is required but not installed.")
  }

  library(foreach)
  library(doParallel)

  # ------------------------------ 0) Setup Clustering -------------------------
  cluster <- if ("c" %in% names(s)) "c" else NULL

  # ------------------------------ 1) Checks -----------------------------------
  stopifnot(is.data.frame(s))

  # Allowed methods (removed "em"):
  valid_methods <- c("ols1", "ols2", "did", "ebal", "aipw")
  if (!(method %in% valid_methods)) {
    stop("method must be one of: ", paste(valid_methods, collapse=", "))
  }

  # Allowed variance types:
  valid_vtypes <- c("robust","bootstrap","jackknife")
  if (!(vartype %in% valid_vtypes)) {
    stop("vartype must be one of: ", paste(valid_vtypes, collapse=", "))
  }

  if (!"G" %in% names(s)) {
    stop("No column named 'G' found in s (the treatment indicator).")
  }

  # Covariates are assumed to be named x1, x2, ...
  covar <- grep("^x\\d+$", names(s), value=TRUE)  # local variable

  # If ref_period has multiple years, pick only the last as baseline
  if (length(ref_period) > 1) {
    ref_period <- max(ref_period)
  }

  # Possibly drop missing rows in relevant columns
  needed_cols <- c("G", covar)
  if (!is.null(cluster)) needed_cols <- c(needed_cols, cluster)
  if (na.rm) {
    s <- s[complete.cases(s[, needed_cols, drop=FALSE]), ]
  }

  # Identify numeric columns that might be year outcomes
  all_cols_numeric <- names(s)[sapply(s, is.numeric)]
  exclude <- c("G", covar)
  if (!is.null(cluster)) exclude <- c(exclude, cluster)
  numeric_year_cols <- setdiff(all_cols_numeric, exclude)

  numeric_years <- suppressWarnings(as.numeric(numeric_year_cols))
  numeric_years <- numeric_years[!is.na(numeric_years)]
  numeric_years <- sort(unique(numeric_years))

  # Check user-supplied year columns
  Y_tr_cols <- as.character(tr_period)
  Y_ref_col <- as.character(ref_period)
  missing_tr <- setdiff(Y_tr_cols, numeric_year_cols)
  if (length(missing_tr) > 0) {
    stop("Some treatment-year columns not found in data: ",
         paste(missing_tr, collapse=", "))
  }
  if (!(Y_ref_col %in% numeric_year_cols)) {
    stop("Reference-year column '", Y_ref_col, "' not found in data.")
  }

  # ----------------- 2) Single-Period FDID (static) --------------------------
  # Create a single outcome: average over tr_period minus ref_period
  s$tempY <- rowMeans(s[, Y_tr_cols, drop=FALSE], na.rm=TRUE) - s[[Y_ref_col]]

  # If na.rm=TRUE, drop missing rows relevant for that outcome
  needed_for_outcome <- c("tempY", needed_cols)
  if (na.rm) {
    s <- s[complete.cases(s[, needed_for_outcome, drop=FALSE]), ]
  }

  # -------------- Subroutines for each method --------------------------------

  # did: Simple difference in means or cluster-robust difference in means
  est_did <- function(d, covar) {
    fml <- tempY ~ G
    if (is.null(cluster)) {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata")
    } else {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata",
                                 clusters=d[[cluster]])
    }
    cfit <- summary(mod)$coefficients["G", c(1,2,5,6)]
    names(cfit) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    cfit
  }

  # ols1: Include covariates linearly; standard robust or cluster-robust SE
  est_ols1 <- function(d, covar) {
    if (length(covar) > 0) {
      fml <- as.formula(paste("tempY ~ G +", paste(covar, collapse=" + ")))
    } else {
      fml <- tempY ~ G
    }
    if (is.null(cluster)) {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata")
    } else {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata",
                                 clusters=d[[cluster]])
    }
    cfit <- summary(mod)$coefficients["G", c(1,2,5,6)]
    names(cfit) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    cfit
  }

  # ols2: Covariates plus "super-population" correction if no cluster
  est_ols2 <- function(d, covar) {
    # If no covariates, fallback to simple DID
    if (length(covar) == 0) {
      return(est_did(d, covar))
    }
    # If covariates exist:
    X <- scale(d[, covar, drop=FALSE])
    Y <- d$tempY
    Z <- d$G
    n <- nrow(X)
    p <- ncol(X)

    if (is.null(cluster)) {
      # Fit full interaction: Y ~ Z*X
      linmod <- lm(Y ~ Z * X)
      est <- coef(linmod)[2]
      vmat <- car::hccm(linmod)  # robust var-cov matrix
      vehw <- vmat[2,2]

      # Interaction terms: (p+3) : (2p+2)
      inter <- coef(linmod)[(p + 3):(2 * p + 2)]
      var_super <- vehw + sum(inter * (cov(X) %*% inter)) / n
      se_super  <- sqrt(var_super)

      out <- c(est,
               se_super,
               est - 1.96 * se_super,
               est + 1.96 * se_super)
    } else {
      df_lin <- data.frame(tempY = Y, Z = Z, X)
      fml <- tempY ~ Z * .
      mod <- estimatr::lm_robust(fml, data = df_lin,
                                 se_type = "stata",
                                 clusters = d[[cluster]])
      out <- summary(mod)$coefficients["Z", c(1, 2, 5, 6)]
    }
    names(out) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    out
  }

  # ebal: Entropy Balancing
  est_ebal <- function(d, covar) {
    if (!requireNamespace("ebal", quietly=TRUE)) {
      stop("Method 'ebal' requires the 'ebal' package.")
    }
    # If no covariates, ebal can't balance anything; just do a simple DID:
    if (length(covar) == 0) {
      return(est_did(d, covar))
    }

    # Prepare Treatment indicator as logical
    Treatment <- d$G == 1

    # Covariates
    X <- d[, covar, drop=FALSE]

    # Apply ebal
    eout <- ebal::ebalance(Treatment = Treatment, X = X)

    # Combine weights: 1 for treated, eout$w for control
    wts <- rep(0, nrow(d))
    wts[Treatment] <- 1
    wts[!Treatment] <- eout$w

    # Weighted regression to get robust or cluster-robust SEs
    fml <- tempY ~ G
    if (is.null(cluster)) {
      mod <- estimatr::lm_robust(fml, data=d, weights=wts, se_type="stata")
    } else {
      mod <- estimatr::lm_robust(fml, data=d, weights=wts, se_type="stata",
                                 clusters=d[[cluster]])
    }
    cfit <- summary(mod)$coefficients["G", c(1,2,5,6)]
    names(cfit) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    cfit
  }

  # AIPW estimation using grf::causal_forest
  est_aipw <- function(d, covar) {
    if (!requireNamespace("grf", quietly=TRUE)) {
      stop("Method 'aipw' requires the 'grf' package.")
    }
    X <- d[, covar, drop=FALSE]
    Y <- d$tempY
    W <- d$G
    c.forest <- grf::causal_forest(X, Y, W, seed=1234)
    att.obj  <- grf::average_treatment_effect(
      c.forest, target.sample="treated", method="AIPW"
    )
    est   <- att.obj[1]
    se_   <- att.obj[2]
    lower <- est - 1.96*se_
    upper <- est + 1.96*se_
    c(Estimate=est, Std.Error=se_, CI_Lower=lower, CI_Upper=upper)
  }

  # Master dispatcher
  run_method <- function(d, method, covar) {
    if (method=="did")  return(est_did(d, covar))
    if (method=="ols1") return(est_ols1(d, covar))
    if (method=="ols2") return(est_ols2(d, covar))
    if (method=="ebal") return(est_ebal(d, covar))
    if (method=="aipw") return(est_aipw(d, covar))
    stop("Unknown method: ", method)
  }

  # ------------------ Variance wrappers --------------------------------------

  # 1) "robust" = run_method(d, method, covar) directly
  do_robust <- function(d, method, covar) {
    run_method(d, method, covar)
  }

  # 2) "bootstrap"
  do_bootstrap <- function(d, method, covar, B = nsims) {
    main_est <- run_method(d, method, covar)["Estimate"]

    one_boot <- function(i) {
      idx <- sample.int(nrow(d), replace = TRUE)
      db  <- d[idx, ]
      run_method(db, method, covar)["Estimate"]
    }

    if (!parallel) {
      vals <- replicate(B, one_boot(1), simplify=TRUE)
    } else {
      if (!requireNamespace("foreach", quietly=TRUE) ||
          !requireNamespace("doParallel", quietly=TRUE)) {
        stop("Parallel requires 'foreach' and 'doParallel' packages installed.")
      }
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      vals <- foreach::foreach(
        i = 1:B,
        .combine = c,
        .export = c("run_method",
                    "est_did", "est_ols1", "est_ols2", "est_ebal", "est_aipw")
      ) %dopar% {
        one_boot(i)
      }
      parallel::stopCluster(cl)
    }
    se_ <- stats::sd(vals, na.rm = TRUE)
    ci_ <- stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)
    c(Estimate = main_est, Std.Error = se_, CI_Lower = ci_[1], CI_Upper = ci_[2])
  }

  # 3) "jackknife"
  do_jackknife <- function(d, method, covar) {
    n <- nrow(d)
    main_est <- run_method(d, method, covar)["Estimate"]

    one_jack <- function(i) {
      d_j <- d[-i, ]
      run_method(d_j, method, covar)["Estimate"]
    }

    if (!parallel) {
      jvals <- numeric(n)
      for (i in seq_len(n)) {
        jvals[i] <- one_jack(i)
      }
    } else {
      if (!requireNamespace("foreach", quietly=TRUE) ||
          !requireNamespace("doParallel", quietly=TRUE)) {
        stop("Parallel requires 'foreach' and 'doParallel' packages installed.")
      }
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      jvals <- foreach::foreach(
        i = 1:n,
        .combine = c,
        .export = c("run_method",
                    "est_did","est_ols1","est_ols2","est_ebal","est_aipw")
      ) %dopar% {
        one_jack(i)
      }
      parallel::stopCluster(cl)
    }

    est_mean <- mean(jvals, na.rm = TRUE)
    # Jackknife SE formula
    se_ <- sqrt((n - 1)/n * sum((jvals - est_mean)^2))
    ci_ <- stats::quantile(jvals, probs = c(0.025, 0.975), na.rm = TRUE)
    c(Estimate = main_est, Std.Error = se_, CI_Lower = ci_[1], CI_Upper = ci_[2])
  }

  get_estimate <- function(d, method, vartype, covar) {
    if (vartype == "robust") {
      do_robust(d, method, covar)
    } else if (vartype == "bootstrap") {
      do_bootstrap(d, method, covar, nsims)
    } else if (vartype == "jackknife") {
      do_jackknife(d, method, covar)
    }
  }

  # Single FDID estimate (static)
  static_result_vec <- get_estimate(s, method, vartype, covar)
  static_result <- as.data.frame(t(static_result_vec))

  # ------------- 3) Dynamic FDID for each year column ------------------------
  dynamic_df <- data.frame(
    Estimate  = numeric(length(numeric_years)),
    Std.Error = numeric(length(numeric_years)),
    CI_Lower  = numeric(length(numeric_years)),
    CI_Upper  = numeric(length(numeric_years)),
    row.names = as.character(numeric_years)
  )
  last_ref_num <- as.numeric(ref_period)

  for (i in seq_along(numeric_years)) {
    t_yr <- numeric_years[i]
    if (t_yr == last_ref_num) {
      # By convention, the effect at the reference year is 0
      dynamic_df[i, ] <- 0
    } else {
      t_col <- as.character(t_yr)
      ref_col <- as.character(last_ref_num)
      d_dyn <- s
      d_dyn$tempY <- d_dyn[[t_col]] - d_dyn[[ref_col]]
      if (na.rm) {
        needed_dyn <- c("tempY", "G")
        if (!is.null(cluster)) needed_dyn <- c(needed_dyn, cluster)
        needed_dyn <- c(needed_dyn, covar)
        d_dyn <- d_dyn[complete.cases(d_dyn[, needed_dyn, drop=FALSE]), ]
      }
      est_vec <- get_estimate(d_dyn, method, vartype, covar)
      dynamic_df[i, ] <- est_vec
    }
  }

  # ------------- 4) Raw Means by group & year --------------------------------
  rawdf <- data.frame(
    year  = numeric(0),
    group = character(0),
    meanY = numeric(0)
  )
  for (yr in numeric_years) {
    col_yr <- as.character(yr)
    if (!col_yr %in% names(s)) {
      warning("Column '", col_yr, "' not found in data. Skipping raw_means.")
      next
    }
    rawdf <- rbind(rawdf,
                   data.frame(
                     year  = yr,
                     group = "Treated",
                     meanY = mean(s[s$G == 1, col_yr][[col_yr]], na.rm=TRUE)
                   ),
                   data.frame(
                     year  = yr,
                     group = "Control",
                     meanY = mean(s[s$G == 0, col_yr][[col_yr]], na.rm=TRUE)
                   ))
  }

  # ------------- 5) Optional: Propensity Scores (for "aipw") -----------------
  ps_vec <- NULL
  if (method == "aipw") {
    if (length(covar) == 0) {
      # If no covariates, just use overall mean as the "propensity"
      ps_vec <- rep(mean(s$G, na.rm=TRUE), nrow(s))
    } else {
      ps_formula <- as.formula(paste("G ~", paste(covar, collapse=" + ")))
      ps_fit <- stats::glm(ps_formula, data=s, family=binomial)
      ps_vec <- stats::predict(ps_fit, newdata=s, type="response")
    }
  }

  # ------------- 6) Compile & Return -----------------------------------------
  out <- list(
    static     = static_result,
    dynamic    = dynamic_df,
    raw_means  = rawdf,
    tr_period  = tr_period,
    ref_period = ref_period,
    method     = method,
    vartype    = vartype,
    years      = numeric_years,
    G          = s$G,
    ps         = ps_vec,
    call       = match.call()
  )
  class(out) <- "fdid"
  return(out)
}
