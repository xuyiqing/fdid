#' Plot Multiple FDID Estimates
#'
#' Creates a side-by-side comparison plot of point estimates and confidence
#' intervals for every element of an `fdid_list`.
#'
#' @param x An object of class `"fdid_list"`.
#' @param palette A palette name from **RColorBrewer**. Defaults to `"Set2"`.
#' @param xlab, main Axis label and title.  If `NULL`, sensible defaults are
#'   used.
#' @param ylim Optional numeric vector of length two giving the y-axis limits.
#' @param ... Additional graphics parameters passed to `plot()`.
#'
#' @return Invisibly returns `x`; called for its side-effect of drawing a plot.
#'
#' @author Rivka Lipkovitz
#'
#' @export
plot.fdid_list <- function(x,
                           palette = "Set2",
                           xlab    = NULL,
                           main    = NULL,
                           ylim    = NULL,
                           ...) {

  if (!inherits(x, "fdid_list")) {
    stop("Input must be of class 'fdid_list'.")
  }
  x <- rev(x)

  # ---- Dependencies & palette checks ---------------------------------------
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      "Package 'RColorBrewer' is required. ",
      "Install it with install.packages('RColorBrewer')."
    )
  }
  pal_info <- RColorBrewer::brewer.pal.info
  if (!palette %in% rownames(pal_info)) {
    stop("Palette '", palette, "' not found. Choose one of: ",
         paste(rownames(pal_info), collapse = ", "))
  }
  n  <- length(x)
  max_cols <- pal_info[palette, "maxcolors"]
  if (n > max_cols) {
    stop("Palette '", palette, "' supports up to ",
         max_cols, " colours; you have ", n, " objects.")
  }
  cols <- RColorBrewer::brewer.pal(max(n, 3), palette)[seq_len(n)]

  # ---- Extract estimates ----------------------------------------------------
  ests <- do.call(
    rbind,
    lapply(x, function(obj) {
      # Expect the object to store method and event-study estimate table
      if (is.null(obj$method) || is.null(obj$est$event)) {
        stop("Each 'fdid' object must contain $method and $est$event elements.")
      }
      ev <- obj$est$event
      if (!all(c("Estimate", "CI_Lower", "CI_Upper") %in% names(ev))) {
        stop("Each $est$event must contain 'Estimate', 'CI_Lower', 'CI_Upper'.")
      }
      data.frame(
        method    = obj$method,
        Estimate  = ev$Estimate,
        CI_Lower  = ev$CI_Lower,
        CI_Upper  = ev$CI_Upper,
        stringsAsFactors = FALSE
      )
    })
  )

  # ---- Axis limits ----------------------------------------------------------
  if (is.null(ylim)) {
    rng  <- range(ests$CI_Lower, ests$CI_Upper, na.rm = TRUE)
    ylim <- rng * 1.05
  }

  # ---- Base plot ------------------------------------------------------------
  old_mar <- par("mar")
  on.exit(par(mar = old_mar), add = TRUE)
  par(mar = c(4, 4, 2, 2))                # give extra room for y labels

  plot(NA,
       xlim = ylim,
       ylim = c(0.5, n + 0.5),
       yaxt = "n",
       xlab = if (is.null(xlab)) "Estimate" else xlab,
       ylab = if (is.null(xlab)) "Method" else xlab,
       main = if (is.null(main)) "Comparison of FDID Estimates" else main,
       ...)

  abline(v = 0, lty = 2)

  axis(2, at = seq_len(n), labels = ests$method, las = 1)

  # ---- Draw CIs and points --------------------------------------------------
  for (i in seq_len(n)) {
    segments(ests$CI_Lower[i], i, ests$CI_Upper[i], i,
             lwd = 2, col = cols[i])
    points(ests$Estimate[i], i, pch = 16, col = cols[i])
  }

  box()
  invisible(x)
}
