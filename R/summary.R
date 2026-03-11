utils::globalVariables("tr_period")

#' Summary Method for FDID Objects
#'
#' @param object An object of class \code{fdid}.
#' @param ... Additional arguments (not used).
#'
#' @return Prints a summary of the \code{fdid} object.
#' @examples
#' \donttest{
#' data(fdid)
#' mortality$uniqueid <- paste(mortality$provid, mortality$countyid, sep = "-")
#' mortality$G <- ifelse(mortality$pczupu >= median(mortality$pczupu, na.rm = TRUE), 1, 0)
#' s <- fdid_prepare(
#'   data = mortality, Y_label = "mortality",
#'   X_labels = c("avggrain", "lnpop"),
#'   G_label = "G", unit_label = "uniqueid", time_label = "year"
#' )
#' result <- fdid(s, tr_period = 1958:1961, ref_period = 1957)
#' summary(result)
#' }
#' @author Rivka Lipkovitz, Enhan Liu
#' @export
summary.fdid <- function(object, ...) {
  if (!inherits(object, "fdid")) {
    stop("Object must be of class 'fdid'.")
  }


  # --- Helper: format one row ---
  fmt_row <- function(label, est, se, ci_lo, ci_hi, lw = 16, nw = 11) {
    sprintf("  %-*s %*s  %*s  [%s, %s]",
            lw, label,
            nw, formatC(est, format = "f", digits = 4),
            nw, formatC(se, format = "f", digits = 4),
            formatC(ci_lo, format = "f", digits = 4),
            formatC(ci_hi, format = "f", digits = 4))
  }

  fmt_header <- function(lw = 16, nw = 11) {
    sprintf("  %-*s %*s  %*s  %s", lw, "", nw, "Estimate", nw, "Std.Error", "   95% CI")
  }

  sep        <- paste(rep("\u2500", 72), collapse = "")
  sep_double <- paste(rep("\u2550", 72), collapse = "")
  dots       <- paste(rep("\u00B7", 68), collapse = "")

  # --- Period labels ---
  pre_times  <- object$entire_period[object$entire_period < min(object$tr_period)]
  post_times <- object$entire_period[object$entire_period > max(object$tr_period)]

  # === Header ===
  cat("\n")
  cat("  Factorial Difference-in-Differences (FDID) Summary\n")
  cat(" ", sep_double, "\n")
  cat("  Method:            ", object$method, "\n")
  cat("  Variance Type:     ", object$vartype, "\n")
  cat("  Reference Period:  ", paste(object$ref_period, collapse = ", "), "\n")
  if (length(pre_times) > 0)
    cat("  Pre-Event Period:  ", paste(pre_times, collapse = ", "), "\n")
  if (length(object$tr_period) > 0)
    cat("  Event Period:      ", paste(object$tr_period, collapse = ", "), "\n")
  if (length(post_times) > 0)
    cat("  Post-Event Period: ", paste(post_times, collapse = ", "), "\n")
  cat(" ", sep_double, "\n\n")

  # === Aggregate Estimates ===
  cat("  Aggregate Estimates\n")
  cat(" ", sep, "\n")
  cat(fmt_header(), "\n")
  cat(" ", sep, "\n")

  # Helper: check if an estimate row is valid (not all NA)
  is_valid <- function(df) {
    !is.null(df) && !all(is.na(df$Estimate))
  }

  pre  <- object$est$pre
  evt  <- object$est$event
  post <- object$est$post

  if (is_valid(pre))
    cat(fmt_row("Pre-Event",  pre$Estimate,  pre$Std.Error,  pre$CI_Lower,  pre$CI_Upper),  "\n")
  if (is_valid(evt) && any(object$tr_period %in% object$entire_period))
    cat(fmt_row("Event",      evt$Estimate,  evt$Std.Error,  evt$CI_Lower,  evt$CI_Upper),  "\n")
  if (is_valid(post))
    cat(fmt_row("Post-Event", post$Estimate, post$Std.Error, post$CI_Lower, post$CI_Upper), "\n")
  cat(" ", sep, "\n\n")

  # === Dynamic Estimates ===
  cat("  Dynamic Estimates\n")
  cat(" ", sep, "\n")
  cat(fmt_header(), "\n")
  cat(" ", sep, "\n")

  dyn <- object$dynamic
  earliest_event <- min(object$tr_period)
  latest_event   <- max(object$tr_period)

  for (i in seq_len(nrow(dyn))) {
    yr   <- rownames(dyn)[i]
    yr_n <- as.numeric(yr)

    if (yr_n == earliest_event && i > 1) cat("  ", dots, "\n")
    if (yr_n == latest_event + 1)        cat("  ", dots, "\n")

    cat(fmt_row(yr, dyn$Estimate[i], dyn$Std.Error[i],
                dyn$CI_Lower[i], dyn$CI_Upper[i]), "\n")
  }

  cat(" ", sep, "\n")

  invisible(object)
}
