#' Create an `fdid_list` Object
#'
#' Bundles multiple `fdid` objects into a single list with class
#' `"fdid_list"` for convenient collective handling.
#'
#' @param ... One or more objects of class `"fdid"`, or a single list of them.
#' @param validate Logical; if `TRUE` (default) verify each element inherits
#'   from `"fdid"`.
#'
#' @return A list with classes `c("fdid_list", "list")`.
#'
#' @author Rivka Lipkovitz
#'
#' @export
fdid_list <- function(..., validate = TRUE) {
  x <- list(...)

  # Allow a single list input
  if (length(x) == 1L && is.list(x[[1]]) && !inherits(x[[1]], "fdid")) {
    x <- x[[1]]
  }

  if (validate) {
    bad <- !vapply(x, inherits, logical(1), "fdid")
    if (any(bad)) {
      stop(
        "All elements must inherit from class 'fdid'. ",
        "Invalid element(s): ", paste(which(bad), collapse = ", ")
      )
    }
  }

  class(x) <- c("fdid_list", "list")
  x
}
