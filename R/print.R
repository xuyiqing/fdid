#' Print Method for FDID Objects
#'
#' @param x An object of class `fdid`.
#' @param ... Additional arguments (not used).
#'
#' @return Prints a brief overview of the `fdid` object
#' @author Rivka Lipkovitz.
#' @export
print.fdid <- function(x, ...) {
  if (!inherits(x, "fdid")) {
    stop("Object must be of class 'fdid'.")
  }

  cat("FDID Object\n")
  cat("--------------------------------------------------\n")
  cat("Method:           ", x$method, "\n")
  cat("Variance Type:    ", x$vartype, "\n")
  cat("Event Period: ", paste(x$tr_period, collapse = ", "), "\n")
  cat("Reference Period: ", paste(x$ref_period, collapse = ", "), "\n\n")

  cat("Event Period Estimate:\n")
  print(x$est$event)

  invisible(x)
}
