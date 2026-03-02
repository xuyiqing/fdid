#' Plot Multiple FDID Estimates
#'
#' Creates a comparison plot of point estimates and confidence intervals for
#' every element of an `fdid_list`.
#'
#' @param x An object of class `"fdid_list"`.
#' @param xlab,ylab,main Axis labels and title. If `NULL`, sensible defaults are used.
#' @param ylim Optional numeric vector of length two giving the *estimate-axis* limits.
#'   (Backward compatible: for horizontal plots this is the x-limit; for vertical plots this is the y-limit.)
#' @param vertical Logical; default is \code{TRUE}.
#' @param show_vartype Logical; include vartype in labels. Default is \code{TRUE}.
#' @param ... Additional graphics parameters passed to \code{plot()}.
#'
#' @return Invisibly returns `x`; called for its side-effect of drawing a plot.
#'
#' @author Rivka Lipkovitz, Enhan Liu
#'
#' @export
plot.fdid_list <- function(x,
                           xlab    = NULL,
                           ylab    = NULL,
                           main    = NULL,
                           ylim    = NULL,
                           vertical = TRUE,
                           show_vartype = TRUE,
                           ...) {

  if (!inherits(x, "fdid_list")) {
    stop("Input must be of class 'fdid_list'.")
  }


  n <- length(x)
  cols <- rep("black", n)

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

      vt <- obj$vartype
      if (is.null(vt)) vt <- NA_character_

      data.frame(
        method    = as.character(obj$method),
        vartype   = as.character(vt),
        Estimate  = ev$Estimate,
        CI_Lower  = ev$CI_Lower,
        CI_Upper  = ev$CI_Upper,
        stringsAsFactors = FALSE
      )
    })
  )

  # ---- Labels (only show vartype for bootstrap/jackknife; no dot) ----
  vt <- ests$vartype
  lab <- ests$method

  show_vt <- !is.na(vt) & (vt%in% c("bootstrap", "jackknife"))
  lab[show_vt] <- paste(ests$method[show_vt], ests$vartype[show_vt])

  # ---- Axis limits ----------------------------------------------------------
  if (is.null(ylim)) {
    rng  <- range(ests$CI_Lower, ests$CI_Upper, na.rm = TRUE)
    ylim <- rng * 1.05
  }

  # ---- Base plot ------------------------------------------------------------
  old_mar <- par("mar")
  on.exit(par(mar = old_mar), add = TRUE)

  if (is.null(main)) main <- "Comparison of FDID Estimates"

  if (!isTRUE(vertical)) {
    # -------------------- HORIZONTAL (estimate on x-axis) --------------------

    ypos <- rev(seq_len(n))

    par(mar = c(5, 9, 3, 2))  # more room for y tick labels + title space

    plot(NA,
         xlim = ylim,
         ylim = c(0.5, n + 0.5),
         yaxt = "n",
         xlab = if (is.null(xlab)) "Coefficient" else xlab,
         ylab = if (is.null(ylab)) "" else ylab,
         main = main,
         ...)

    abline(v = 0, lty = 2)
    axis(2, at = ypos, labels = lab, las = 1)

    # axis titles moved away to avoid overlap
    #mtext(xlab, side = 1, line = 3)
    #mtext(ylab, side = 2, line = 7)

    for (i in seq_len(n)) {
      segments(ests$CI_Lower[i], ypos[i], ests$CI_Upper[i], ypos[i], lwd = 2, col = cols[i])
      points(ests$Estimate[i], ypos[i], pch = 16, col = cols[i])
    }

  } else {
    # -------------------- VERTICAL (estimates on y-axis) ----------------------

    xpos <- seq_len(n)

    par(mar = c(10, 5, 3, 2))  # big bottom margin for rotated labels

    plot(NA,
         xlim = c(0.5, n + 0.5),
         ylim = ylim,
         xaxt = "n",
         xlab = if (is.null(xlab)) "" else xlab,
         ylab = if (is.null(ylab)) "Coefficient" else ylab,
         main = main,
         ...)

    abline(h = 0, lty = 2)

    axis(1, at = xpos, labels = lab, las = 2)

    # axis titles moved away to avoid overlap
    #mtext(xlab, side = 1, line = 8)
    #mtext(ylab, side = 2, line = 3)

    for (i in seq_len(n)) {
      segments(xpos[i], ests$CI_Lower[i], xpos[i], ests$CI_Upper[i], lwd = 2, col = cols[i])
      points(xpos[i], ests$Estimate[i], pch = 16, col = cols[i])
    }
  }

  box()
  invisible(x)
}
