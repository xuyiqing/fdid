#' Summary Method for FDID Objects
#'
#' @param object An object of class `fdid`.
#' @param ... Additional arguments (not used).
#'
#' @return Prints a summary of the `fdid` object.
#' @author Rivka Lipkovitz
#' @export
summary.fdid <- function(object, ...) {
  if (!inherits(object, "fdid")) {
    stop("Object must be of class 'fdid'.")
  }

  cat("Factorial Difference-in-Differences (FDID) Summary\n")
  cat("--------------------------------------------------\n")
  cat("Method:           ", object$method, "\n")
  cat("Variance Type:    ", object$vartype, "\n")
  cat("Reference Period: ", paste(object$ref_period, collapse = ", "), "\n\n")
  cat("Pre-Event Period: ", paste(object$entire_period[object$entire_period < object$tr_period[1]], collapse = ", "), "\n")
  cat("Event Period: ", paste(object$tr_period, collapse = ", "), "\n")
  cat("Post-Event Period: ", paste(object$entire_period[object$entire_period > object$tr_period[length(tr_period)]], collapse = ", "), "\n")
  cat("Pre-Event Period Estimate:\n")
  print(object$est$pre)
  cat("Event Period Estimate:\n")
  print(object$est$event)
  cat("Post-Event Period Estimate:\n")
  print(object$est$post)
  cat("\nDynamic Estimates:\n")
  print(object$dynamic)

  invisible(object)
}
