#' Factorial Difference-in-Differences Estimation
#'
#' Performs factorial difference-in-differences (FDID) estimation using various methods and variance estimation techniques.
#'
#' @param s A data frame prepared using `fdid_prepare`.
#' @param tr_period A numeric vector specifying the treatment periods.
#' @param ref_period A numeric scalar specifying the reference period.
#' @param entire_period A numeric vector specifying the total range of time periods.
#'   If `NULL`, estimation is performed on all available time periods.
#'   Example: `c(1958, 1959, 1960, 1961)`.
#' @param method A string specifying the estimation method.
#'   Options: `"ols1"`, `"ols2"`, `"did"`, `"ebal"`, `"aipw"`.
#'   Default is `"ols1"`.
#' @param vartype A string specifying the variance estimation type.
#'   Options: `"robust"`, `"bootstrap"`, `"jackknife"`.
#'   Default is `"robust"`.
#' @param missing_data How to handle missing data. Two options:
#'   - `"listwise"`: Drop any row missing **any** relevant column (including
#'     outcomes in the periods used).
#'   - `"available"`: Drop rows only if they are missing in group/covariates/cluster
#'     columns, but allow partial usage of outcomes.
#'   Default is `"listwise"`.
#' @param nsims Number of simulations for bootstrap variance estimation.
#'   Default is `1000`.
#' @param parallel Logical; whether to perform parallel computations.
#'   Default is `FALSE`.
#' @param cores Number of cores for parallel computations.
#'   Default is `2`.
#' @param ATT Logical; if `TRUE`, attempts to estimate the ATT by targeting
#'   the covariate distribution of the treated group (`G=1`). If `FALSE`,
#'   estimates the ATE. Not all methods can handle both; see Details.
#'   Default is `FALSE`.
#'
#' @return A list with the following components:
#'   \item{est}{A list with three elements:
#'              \code{$pre}, \code{$event}, and \code{$post} containing
#'              aggregated pre-treatment, overall event, and post-treatment
#'              FDID estimates, respectively.}
#'   \item{dynamic}{Dynamic FDID estimates for each time in `entire_period`.}
#'   \item{raw_means}{Raw mean outcomes by group for each time in `entire_period`.}
#'   \item{tr_period}{Treatment periods used.}
#'   \item{ref_period}{Reference period used.}
#'   \item{entire_period}{All time periods for dynamic estimation.}
#'   \item{method}{Method used.}
#'   \item{vartype}{Variance type used.}
#'   \item{times}{All numeric time columns found.}
#'   \item{G}{Group indicator (0/1).}
#'   \item{ps}{Propensity scores (if `aipw` method used).}
#'   \item{call}{The matched call.}
#'   \item{ATT}{Logical indicating whether ATT was requested.}
#'
#' @import dplyr
#' @import estimatr
#' @importFrom grf causal_forest
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
#' @author Rivka Lipkovitz
#' @export
fdid <- function(s,
                 tr_period,
                 ref_period,
                 entire_period = NULL,
                 method   = "ols1",
                 vartype  = "robust",
                 missing_data = c("listwise", "available"),
                 nsims    = 1000,
                 parallel = FALSE,
                 cores    = 2,
                 ATT      = FALSE) {

  # Helper to suppress ebal output:
  silent_ebalance <- function(Treatment, X) {
    tf_con <- file(tempfile(), open = "wt")
    on.exit({
      sink(type = "message"); sink(type = "output")
      close(tf_con)
    }, add = TRUE)
    sink(tf_con, type = "output"); sink(tf_con, type = "message")
    old_warn <- getOption("warn"); options(warn = -1)

    out <- ebal::ebalance(Treatment = Treatment, X = X)

    options(warn = old_warn)
    return(out)
  }

  # Match arg:
  missing_data <- match.arg(missing_data)

  # --- Basic checks ---
  if (!requireNamespace("estimatr", quietly=TRUE)) {
    stop("Package 'estimatr' not installed. Please install it.")
  }
  if (!requireNamespace("car", quietly=TRUE)) {
    stop("Package 'car' not installed. Please install it.")
  }
  library(foreach)
  library(doParallel)

  # If there's a cluster column, store its name
  cluster <- if ("c" %in% names(s)) "c" else NULL

  stopifnot(is.data.frame(s))

  valid_methods <- c("ols1", "ols2", "did", "ebal", "aipw")
  if (!(method %in% valid_methods)) {
    stop("method must be one of: ", paste(valid_methods, collapse=", "))
  }
  if (method %in% c("did","ols1") && ATT) {
    stop("ATT = TRUE is not supported for method = 'did'.")
  }

  valid_vtypes <- c("robust","bootstrap","jackknife")
  if (!(vartype %in% valid_vtypes)) {
    stop("vartype must be one of: ", paste(valid_vtypes, collapse=", "))
  }

  if (!"G" %in% names(s)) {
    stop("No column named 'G' found in s (the group indicator).")
  }

  # Covariates are those named x1, x2, ... (per your original code)
  covar <- grep("^x\\d+$", names(s), value=TRUE)

  # If user gave multiple ref_period, pick last
  if (length(ref_period) > 1) {
    ref_period <- max(ref_period)
  }

  # Identify numeric columns that might be outcomes
  all_cols_numeric <- names(s)[sapply(s, is.numeric)]
  exclude <- c("G", covar)
  if (!is.null(cluster)) exclude <- c(exclude, cluster)
  numeric_time_cols <- setdiff(all_cols_numeric, exclude)
  numeric_times <- suppressWarnings(as.numeric(numeric_time_cols))
  numeric_times <- numeric_times[!is.na(numeric_times)]
  numeric_times <- sort(unique(numeric_times))

  # Check user-supplied time columns
  Y_tr_cols <- as.character(tr_period)
  Y_ref_col <- as.character(ref_period)
  missing_tr <- setdiff(Y_tr_cols, numeric_time_cols)
  if (length(missing_tr) > 0) {
    stop("Some treatment-time columns not found in data: ",
         paste(missing_tr, collapse=", "))
  }
  if (!(Y_ref_col %in% numeric_time_cols)) {
    stop("Reference-time column '", Y_ref_col, "' not found in data.")
  }

  # Missing-data approach
  needed_cols_minimal <- c("G", covar)
  if (!is.null(cluster)) needed_cols_minimal <- c(needed_cols_minimal, cluster)

  # If entire_period is NULL, use all numeric_times
  if (!is.null(entire_period)) {
    all_times <- sort(intersect(numeric_times, entire_period))
  } else {
    all_times <- numeric_times
  }

  if (missing_data == "listwise") {
    needed_cols_all <- c(needed_cols_minimal, as.character(all_times))
    s <- s[complete.cases(s[, needed_cols_all, drop=FALSE]), ]
  } else {
    # "available": only remove missing in G, covars, cluster
    s <- s[complete.cases(s[, needed_cols_minimal, drop=FALSE]), ]
  }

  # --- AUTOMATICALLY ENCODE ANY FACTOR COVARIATES FOR AIPW ---
  if (method == "aipw" && length(covar) > 0) {
    # Check if any of these covars are not numeric
    is_non_num <- sapply(s[, covar, drop=FALSE], function(z) !is.numeric(z))
    if (any(is_non_num)) {
      # Build a formula for model.matrix
      form_covar <- as.formula(paste("~", paste(covar, collapse = " + ")))
      mm <- model.matrix(form_covar, data = s)
      mm <- mm[, -1, drop = FALSE]  # drop the intercept column

      # Drop old covar columns
      keep_cols <- setdiff(names(s), covar)
      # Combine
      s <- cbind(s[, keep_cols, drop=FALSE], mm)
      # Update covar to match new dummy columns
      covar <- colnames(mm)
    }
  }

  # Create single FDID outcome: average of treat periods minus reference
  s$tempY <- rowMeans(s[, Y_tr_cols, drop=FALSE], na.rm=TRUE) - s[[Y_ref_col]]

  # ----- Estimation Subroutines -----

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

  est_ols1 <- function(d, covar) {
    # Y ~ G + X (no interactions)
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

  est_ols2 <- function(d, covar) {
    # Y ~ G * X (with interactions).
    if (length(covar) == 0) {
      return(est_did(d, covar))
    }

    # Prep data
    X <- scale(d[, covar, drop = FALSE])
    Y <- d$tempY
    Z <- d$G
    n <- nrow(X)

    # Fit a plain lm to get coefficients
    linmod <- lm(Y ~ Z * X)
    coefs  <- coef(linmod)

    # Names of the interaction terms ("Z:X1", "Z:X2", …)
    int_names <- paste0("Z:X", colnames(X))

    # Mean of each (scaled) covariate in the treated group
    X_bar <- colMeans(X[Z == 1, , drop = FALSE])

    if (ATT) {
      # — Point estimate at X̄ :
      beta_Z  <- coefs["Z"]
      beta_ZX <- coefs[int_names]
      est_super <- beta_Z + sum(beta_ZX * X_bar)

      # — Variance matrix
      if (is.null(cluster)) {
        V <- car::hccm(linmod)                      # HC0
      } else {
        # require sandwich for cluster‐robust
        if (!requireNamespace("sandwich", quietly = TRUE)) {
          stop("Please install the 'sandwich' package for clustered SEs")
        }
        V <- sandwich::vcovCL(linmod, cluster = d[[cluster]])
      }

      all_nms <- names(coefs)
      g <- numeric(length(coefs)); names(g) <- all_nms
      g["Z"] <- 1
      g[int_names] <- X_bar

      var_super <- as.numeric(t(g) %*% V %*% g)
      se_super  <- sqrt(var_super)

      out <- c(
        Estimate  = est_super,
        Std.Error = se_super,
        CI_Lower  = est_super - 1.96 * se_super,
        CI_Upper  = est_super + 1.96 * se_super
      )

    } else {
      # — classic effect at X=0
      est0 <- coefs["Z"]

      if (is.null(cluster)) {
        V0  <- car::hccm(linmod)
        se0 <- sqrt(V0["Z","Z"])
      } else {
        if (!requireNamespace("sandwich", quietly = TRUE)) {
          stop("Please install the 'sandwich' package for clustered SEs")
        }
        V0  <- sandwich::vcovCL(linmod, cluster = d[[cluster]])
        se0 <- sqrt(V0["Z","Z"])
      }

      out <- c(
        Estimate  = est0,
        Std.Error = se0,
        CI_Lower  = est0 - 1.96 * se0,
        CI_Upper  = est0 + 1.96 * se0
      )
    }

    names(out) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    out
  }

  est_ebal <- function(d, covar) {
    if (!requireNamespace("ebal", quietly=TRUE)) {
      stop("Method 'ebal' requires the 'ebal' package.")
    }
    # ebal code that reweights controls to match the distribution of treated => ATT
    if (!ATT) {
      stop("Currently, 'ebal' only implemented for ATT=TRUE. Please set ATT=TRUE or choose another method.")
    }
    if (length(covar) == 0) {
      return(est_did(d, covar))
    }
    group1 <- d$G == 1
    X <- d[, covar, drop=FALSE]
    eout <- silent_ebalance(Treatment = group1, X = X)

    wts <- rep(0, nrow(d))
    wts[group1] <- 1
    wts[!group1] <- eout$w

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
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Method 'aipw' requires the 'grf' package.")
    }
    X <- d[, covar, drop = FALSE]
    Y <- d$tempY
    W <- d$G

    c.forest <- grf::causal_forest(X, Y, W)
    att      <- grf::average_treatment_effect(
      c.forest,
      target.sample = ifelse(ATT, "treated", "all"),
      method        = "AIPW"
    )
    est <- as.numeric(att["estimate"])
    se_ <- as.numeric(att["std.err"])
    c(
      Estimate  = est,
      Std.Error = se_,
      CI_Lower  = est - 1.96 * se_,
      CI_Upper  = est + 1.96 * se_
    )
  }


  run_method <- function(d, method, covar) {
    if (method=="did")  return(est_did(d, covar))
    if (method=="ols1") return(est_ols1(d, covar))
    if (method=="ols2") return(est_ols2(d, covar))
    if (method=="ebal") return(est_ebal(d, covar))
    if (method=="aipw") return(est_aipw(d, covar))
    stop("Unknown method: ", method)
  }

  # ----- Variance wrappers -----

  do_robust <- function(d, method, covar) {
    run_method(d, method, covar)
  }

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
        stop("Parallel requires 'foreach' and 'doParallel' packages.")
      }
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      vals <- foreach::foreach(
        i = 1:B,
        .combine = c,
        .export = c("run_method","est_did", "est_ols1", "est_ols2", "est_ebal", "est_aipw")
      ) %dopar% {
        one_boot(i)
      }
      parallel::stopCluster(cl)
    }
    se_ <- sd(vals, na.rm=TRUE)
    ci_ <- quantile(vals, probs=c(0.025, 0.975), na.rm=TRUE)
    c(Estimate=main_est, Std.Error=se_, CI_Lower=ci_[1], CI_Upper=ci_[2])
  }

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
        stop("Parallel requires 'foreach' and 'doParallel' packages.")
      }
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      jvals <- foreach::foreach(
        i = 1:n,
        .combine = c,
        .export = c("run_method","est_did","est_ols1","est_ols2","est_ebal","est_aipw")
      ) %dopar% {
        one_jack(i)
      }
      parallel::stopCluster(cl)
    }
    est_mean <- mean(jvals, na.rm=TRUE)
    se_ <- sqrt((n - 1)/n * sum((jvals - est_mean)^2))
    ci_ <- quantile(jvals, probs=c(0.025, 0.975), na.rm=TRUE)
    c(Estimate=main_est, Std.Error=se_, CI_Lower=ci_[1], CI_Upper=ci_[2])
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

  # ----- Single FDID estimate (static) -----
  static_result_vec <- get_estimate(s, method, vartype, covar)
  static_result <- as.data.frame(t(static_result_vec))

  # ----- Dynamic FDID for each time column -----
  dynamic_times <- all_times
  dynamic_df <- data.frame(
    Estimate  = numeric(length(dynamic_times)),
    Std.Error = numeric(length(dynamic_times)),
    CI_Lower  = numeric(length(dynamic_times)),
    CI_Upper  = numeric(length(dynamic_times)),
    row.names = as.character(dynamic_times)
  )
  last_ref_num <- as.numeric(ref_period)

  for (i in seq_along(dynamic_times)) {
    t_yr <- dynamic_times[i]
    if (t_yr == last_ref_num) {
      # By definition, the reference period effect is 0
      dynamic_df[i, ] <- 0
    } else {
      t_col <- as.character(t_yr)
      ref_col <- as.character(last_ref_num)
      d_dyn <- s
      d_dyn$tempY <- d_dyn[[t_col]] - d_dyn[[ref_col]]
      d_dyn <- d_dyn[!is.na(d_dyn$tempY), , drop=FALSE]
      est_vec <- get_estimate(d_dyn, method, vartype, covar)
      dynamic_df[i, ] <- est_vec
    }
  }

  # ----- Raw means by group & time -----
  rawdf <- data.frame(
    time  = numeric(0),
    group = character(0),
    meanY = numeric(0)
  )
  for (yr in dynamic_times) {
    col_yr <- as.character(yr)
    if (!col_yr %in% names(s)) {
      warning("Column '", col_yr, "' not found in data. Skipping raw_means.")
      next
    }
    rawdf <- rbind(rawdf,
                   data.frame(
                     time  = yr,
                     group = "Group 1",
                     meanY = mean(s[[col_yr]][ s$G == 1 ], na.rm = TRUE)
                   ),
                   data.frame(
                     time  = yr,
                     group = "Group 0",
                     meanY = mean(s[[col_yr]][ s$G == 0 ], na.rm = TRUE)
                   ))
  }

  # ----- Propensity Scores (only for AIPW) -----
  ps_vec <- NULL
  if (method == "aipw") {
    if (length(covar) == 0) {
      ps_vec <- rep(mean(s$G, na.rm=TRUE), nrow(s))
    } else {
      ps_formula <- as.formula(paste("G ~", paste(covar, collapse=" + ")))
      ps_fit <- glm(ps_formula, data=s, family=binomial)
      ps_vec <- predict(ps_fit, newdata=s, type="response")
    }
  }

  # ----- Pre-event & Post-event aggregates -----
  earliest_event_time <- min(tr_period)
  latest_event_time   <- max(tr_period)

  get_aggregate_estimate <- function(sdata, time_vec, ref_col) {
    if (length(time_vec) == 0) {
      return(c(Estimate=NA, Std.Error=NA, CI_Lower=NA, CI_Upper=NA))
    }
    d_agg <- sdata
    d_agg$tempY <- rowMeans(d_agg[, as.character(time_vec), drop=FALSE], na.rm=TRUE) -
      d_agg[[ref_col]]
    d_agg <- d_agg[!is.na(d_agg$tempY), , drop=FALSE]
    get_estimate(d_agg, method, vartype, covar)
  }

  pre_times <- dynamic_times[dynamic_times < earliest_event_time]
  post_times <- dynamic_times[dynamic_times > latest_event_time]
  pre_event_vec  <- get_aggregate_estimate(s, pre_times, Y_ref_col)
  post_event_vec <- get_aggregate_estimate(s, post_times, Y_ref_col)

  pre_event_result  <- as.data.frame(t(pre_event_vec))
  event_result      <- static_result
  post_event_result <- as.data.frame(t(post_event_vec))

  forced_cols <- c("Estimate", "Std.Error", "CI_Lower", "CI_Upper")
  colnames(pre_event_result)   <- forced_cols
  colnames(event_result)       <- forced_cols
  colnames(post_event_result)  <- forced_cols

  est_list <- list(
    pre   = pre_event_result,
    event = event_result,
    post  = post_event_result
  )

  # Compile output
  out <- list(
    est           = est_list,
    dynamic       = dynamic_df,
    raw_means     = rawdf,
    tr_period     = tr_period,
    ref_period    = ref_period,
    entire_period = dynamic_times,
    method        = method,
    vartype       = vartype,
    times         = numeric_times,
    G             = s$G,
    ps            = ps_vec,
    call          = match.call(),
    ATT           = ATT
  )
  class(out) <- "fdid"
  return(out)
}
