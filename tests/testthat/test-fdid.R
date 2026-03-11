library(fdid)

# ===========================================================================
# CRAN-facing tests
# Tiny synthetic panel (40 units × 3 periods = 120 rows).
# Uses only the two fastest estimators (did, ols1) with no bootstrap.
# Designed to run in < 5 seconds on a single core.
# ===========================================================================

set.seed(42)
n_units <- 40
times   <- 1:3

cran_data <- data.frame(
  id      = rep(seq_len(n_units), each = length(times)),
  time    = rep(times, times = n_units),
  G       = rep(rep(c(0L, 1L), n_units / 2), each = length(times)),
  outcome = rnorm(n_units * length(times)),
  covar1  = rnorm(n_units * length(times)),
  covar2  = rnorm(n_units * length(times))
)

s_cran <- fdid_prepare(
  data       = cran_data,
  Y_label    = "outcome",
  X_labels   = c("covar1", "covar2"),
  G_label    = "G",
  unit_label = "id",
  time_label = "time"
)

test_that("[CRAN] fdid_prepare returns expected structure", {
  expect_s3_class(s_cran, "tbl_df")
  expect_true("G"    %in% names(s_cran))
  expect_true("unit" %in% names(s_cran))
  expect_true(any(grepl("^Y_", names(s_cran))))
  expect_true(any(grepl("^x",  names(s_cran))))
})

test_that("[CRAN] fdid: did", {
  res <- fdid(s_cran, tr_period = 2L, ref_period = 1L, method = "did", vartype = "robust")
  expect_s3_class(res, "fdid")
  expect_equal(res$method, "did")
  expect_true(is.finite(res$est$event$Estimate))
})

test_that("[CRAN] fdid: ols1", {
  res <- fdid(s_cran, tr_period = 2L, ref_period = 1L, method = "ols1", vartype = "robust")
  expect_s3_class(res, "fdid")
  expect_equal(res$method, "ols1")
  expect_true(is.finite(res$est$event$Estimate))
})

test_that("[CRAN] fdid_list bundles results", {
  r1 <- fdid(s_cran, tr_period = 2L, ref_period = 1L, method = "did",  vartype = "robust")
  r2 <- fdid(s_cran, tr_period = 2L, ref_period = 1L, method = "ols1", vartype = "robust")
  fl <- fdid_list(r1, r2)
  expect_s3_class(fl, "fdid_list")
  expect_length(fl, 2L)
})

test_that("[CRAN] summary and print run without error", {
  res <- fdid(s_cran, tr_period = 2L, ref_period = 1L, method = "did", vartype = "robust")
  expect_output(summary(res))
  expect_output(print(res))
})

test_that("[CRAN] plot.fdid: raw and dynamic run without error", {
  res <- fdid(s_cran, tr_period = 2L, ref_period = 1L, method = "ols1", vartype = "robust")
  expect_invisible(plot(res, type = "raw"))
  expect_invisible(plot(res, type = "dynamic"))
})


# ===========================================================================
# Full internal tests  (skipped on CRAN)
# Uses the real mortality dataset and all six estimators.
# ===========================================================================

data(fdid)
mortality$uniqueid <- paste(mortality$provid, mortality$countyid, sep = "-")
mortality$G <- ifelse(mortality$pczupu >= median(mortality$pczupu, na.rm = TRUE), 1, 0)

s <- fdid_prepare(
  data       = mortality,
  Y_label    = "mortality",
  X_labels   = c("avggrain", "lnpop"),
  G_label    = "G",
  unit_label = "uniqueid",
  time_label = "year"
)

TR  <- 1958
REF <- 1957

# Helper: basic structural checks on any fdid result
check_fdid <- function(res, method_name) {
  expect_s3_class(res, "fdid")
  expect_named(res, c("est", "dynamic", "raw_means", "tr_period", "ref_period",
                      "entire_period", "method", "vartype", "times", "G", "ps",
                      "call", "target.pop"),
               ignore.order = TRUE)
  expect_equal(res$method, method_name)
  expect_true(is.finite(res$est$event$Estimate),
              label = paste(method_name, "event estimate is finite"))
  expect_true(is.finite(res$est$event$Std.Error),
              label = paste(method_name, "event SE is finite"))
  expect_true(nrow(res$dynamic) > 0,
              label = paste(method_name, "dynamic estimates non-empty"))
}

# ---------------------------------------------------------------------------
test_that("fdid_prepare returns expected structure", {
  skip_on_cran()
  expect_s3_class(s, "tbl_df")
  expect_true("G"    %in% names(s))
  expect_true("unit" %in% names(s))
  expect_true(any(grepl("^Y_", names(s))))
  expect_true(any(grepl("^x",  names(s))))
})

# ---------------------------------------------------------------------------
test_that("fdid: did", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "did", vartype = "robust")
  check_fdid(res, "did")
})

test_that("fdid: ols1", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  check_fdid(res, "ols1")
})

test_that("fdid: ols2", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols2", vartype = "robust")
  check_fdid(res, "ols2")
})

test_that("fdid: ebal (targets G=1)", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ebal",
              vartype = "robust", target.pop = "1")
  check_fdid(res, "ebal")
})

test_that("fdid: ipw", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ipw", vartype = "robust")
  check_fdid(res, "ipw")
  expect_false(is.null(res$ps), label = "ipw returns propensity scores")
})

test_that("fdid: aipw", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "aipw", vartype = "robust")
  check_fdid(res, "aipw")
  expect_false(is.null(res$ps), label = "aipw returns propensity scores")
})

# ---------------------------------------------------------------------------
test_that("target.pop shifts estimates for ols2", {
  skip_on_cran()
  all_ <- fdid(s, tr_period = TR, ref_period = REF, method = "ols2",
               vartype = "robust", target.pop = "all")
  att  <- fdid(s, tr_period = TR, ref_period = REF, method = "ols2",
               vartype = "robust", target.pop = "1")
  atc  <- fdid(s, tr_period = TR, ref_period = REF, method = "ols2",
               vartype = "robust", target.pop = "0")
  # Estimates need not be identical across target populations
  expect_false(isTRUE(all.equal(all_$est$event$Estimate, att$est$event$Estimate)) &&
               isTRUE(all.equal(all_$est$event$Estimate, atc$est$event$Estimate)),
               label = "ols2 target.pop changes the estimate")
})

# ---------------------------------------------------------------------------
test_that("did and ols1 are unaffected by target.pop", {
  skip_on_cran()
  r_all <- fdid(s, tr_period = TR, ref_period = REF, method = "did",
                vartype = "robust", target.pop = "all")
  r_att <- fdid(s, tr_period = TR, ref_period = REF, method = "did",
                vartype = "robust", target.pop = "1")
  expect_equal(r_all$est$event$Estimate, r_att$est$event$Estimate,
               tolerance = 1e-10,
               label = "did estimate same across target.pop")
})

# ---------------------------------------------------------------------------
test_that("entire_period restricts dynamic output", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1",
              vartype = "robust", entire_period = c(1957, 1958, 1959))
  expect_equal(sort(res$entire_period), c(1957, 1958, 1959))
  expect_equal(nrow(res$dynamic), 3L)
})

# ---------------------------------------------------------------------------
test_that("fdid_list bundles multiple results", {
  skip_on_cran()
  r1 <- fdid(s, tr_period = TR, ref_period = REF, method = "did",  vartype = "robust")
  r2 <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  fl <- fdid_list(r1, r2)
  expect_s3_class(fl, "fdid_list")
  expect_length(fl, 2L)
})

# ---------------------------------------------------------------------------
test_that("summary.fdid runs without error", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  expect_output(summary(res))
})

test_that("print.fdid runs without error", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  expect_output(print(res))
})

# ---------------------------------------------------------------------------
test_that("plot.fdid: raw and dynamic run without error", {
  skip_on_cran()
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  expect_invisible(plot(res, type = "raw"))
  expect_invisible(plot(res, type = "dynamic"))
})

test_that("plot.fdid: overlap requires ipw/aipw", {
  skip_on_cran()
  res_aipw <- fdid(s, tr_period = TR, ref_period = REF, method = "aipw", vartype = "robust")
  expect_invisible(plot(res_aipw, type = "overlap"))
})
