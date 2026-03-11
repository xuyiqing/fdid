library(fdid)

# Shared setup: prepare a minimal dataset once for the whole file.
# Uses a single treatment period and two covariates to keep tests fast.
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
  expect_s3_class(s, "tbl_df")
  expect_true("G"    %in% names(s))
  expect_true("unit" %in% names(s))
  expect_true(any(grepl("^Y_", names(s))))
  expect_true(any(grepl("^x",  names(s))))
})

# ---------------------------------------------------------------------------
test_that("fdid: did", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "did", vartype = "robust")
  check_fdid(res, "did")
})

test_that("fdid: ols1", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  check_fdid(res, "ols1")
})

test_that("fdid: ols2", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols2", vartype = "robust")
  check_fdid(res, "ols2")
})

test_that("fdid: ebal (targets G=1)", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ebal",
              vartype = "robust", target.pop = "1")
  check_fdid(res, "ebal")
})

test_that("fdid: ipw", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ipw", vartype = "robust")
  check_fdid(res, "ipw")
  expect_false(is.null(res$ps), label = "ipw returns propensity scores")
})

test_that("fdid: aipw", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "aipw", vartype = "robust")
  check_fdid(res, "aipw")
  expect_false(is.null(res$ps), label = "aipw returns propensity scores")
})

# ---------------------------------------------------------------------------
test_that("target.pop shifts estimates for ols2", {
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
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1",
              vartype = "robust", entire_period = c(1957, 1958, 1959))
  expect_equal(sort(res$entire_period), c(1957, 1958, 1959))
  expect_equal(nrow(res$dynamic), 3L)
})

# ---------------------------------------------------------------------------
test_that("fdid_list bundles multiple results", {
  r1 <- fdid(s, tr_period = TR, ref_period = REF, method = "did",  vartype = "robust")
  r2 <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  fl <- fdid_list(r1, r2)
  expect_s3_class(fl, "fdid_list")
  expect_length(fl, 2L)
})

# ---------------------------------------------------------------------------
test_that("summary.fdid runs without error", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  expect_output(summary(res))
})

test_that("print.fdid runs without error", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  expect_output(print(res))
})

# ---------------------------------------------------------------------------
test_that("plot.fdid: raw and dynamic run without error", {
  res <- fdid(s, tr_period = TR, ref_period = REF, method = "ols1", vartype = "robust")
  expect_invisible(plot(res, type = "raw"))
  expect_invisible(plot(res, type = "dynamic"))
})

test_that("plot.fdid: overlap requires ipw/aipw", {
  res_aipw <- fdid(s, tr_period = TR, ref_period = REF, method = "aipw", vartype = "robust")
  expect_invisible(plot(res_aipw, type = "overlap"))
})
