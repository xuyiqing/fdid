#' Plot Results from FDID Analysis
#'
#' Provides visualisations for FDID results, including raw means, dynamic
#' effects, and propensity-score overlap.  The comparison plot of multiple
#' methods has been removed; use `plot.fdid_list()` for that.
#'
#' @param x   An `fdid` object.
#' @param type One of `"raw"`, `"dynamic"`, or `"overlap"`.
#' @param connected Logical; if `TRUE`, connects points with lines in
#'        the "raw" and "dynamic" plots.  Default is `FALSE`.
#' @param alpha_shade Transparency for shading the treatment period.
#' @param palette A palette name from **RColorBrewer**.  Default `"Set2"`.
#' @param group_labels Labels for the two groups.
#' @param xlab,ylab,main Axis labels and main title.
#' @param ylim Y-axis limits.  Default `NULL` (computed automatically).
#' @param ...  Additional graphics parameters.
#'
#' @return Produces a plot; invisibly returns `NULL`.
#' @author Rivka Lipkovitz
#' @export
plot.fdid <- function(x,
                      type = c("raw", "dynamic", "overlap"),
                      connected = FALSE,
                      alpha_shade = 0.2,
                      palette = "Set2",
                      group_labels = c("Group 0", "Group 1"),
                      xlab = NULL,
                      ylab = NULL,
                      main = NULL,
                      ylim = NULL,
                      ...) {

  type <- match.arg(type)

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please install the 'RColorBrewer' package.")
  }

  available_palettes <- rownames(RColorBrewer::brewer.pal.info)
  if (!(palette %in% available_palettes)) {
    stop("Palette not found. Choose from: ",
         paste(available_palettes, collapse = ", "))
  }

  group_colors <- RColorBrewer::brewer.pal(3, palette)[1:2]

  shade_treatment <- function(tr_period, yrange, col = "gray") {
    rect(min(tr_period) - 0.5,
         yrange[1] - 10,
         max(tr_period) + 0.5,
         yrange[2] + 10,
         col = adjustcolor(col, alpha.f = alpha_shade),
         border = NA)
  }

  ## ------------------------------------------------------------------
  ## RAW means
  ## ------------------------------------------------------------------
  if (type == "raw") {
    rawdf <- x$raw_means
    if (is.null(rawdf)) stop("No 'raw_means' found in the 'fdid' object.")
    if (!all(c("time", "group", "meanY") %in% names(rawdf))) {
      stop("'raw_means' is missing required columns.")
    }

    times <- sort(unique(rawdf$time))
    if (is.null(ylim)) {
      tmp <- range(rawdf$meanY, na.rm = TRUE)
      ylim <- c(min(tmp[1], 0.99 * tmp[1]), max(tmp[2], 1.01 * tmp[2]))
    }

    plot(NULL, xlim = range(times), ylim = ylim,
         xlab = xlab %||% "Time",
         ylab = ylab %||% "Mean outcome",
         main = main, mgp = c(2, .8, 0), ...)
    shade_treatment(x$tr_period, ylim)

    g1 <- rawdf[rawdf$group == "Group 1", ]
    g0 <- rawdf[rawdf$group == "Group 0", ]

    if (connected) {
      lines(g1$time, g1$meanY, col = group_colors[2], lwd = 2)
      lines(g0$time, g0$meanY, col = group_colors[1], lwd = 2)
    }
    points(g1$time, g1$meanY, pch = 16, col = group_colors[2])
    points(g0$time, g0$meanY, pch = 16, col = group_colors[1])

    legend("topleft", legend = group_labels,
           col = group_colors, pch = 16, bty = "n")
    return(invisible(NULL))
  }

  ## ------------------------------------------------------------------
  ## DYNAMIC effects
  ## ------------------------------------------------------------------
  if (type == "dynamic") {
    dyn <- x$dynamic
    if (is.null(dyn) || nrow(dyn) == 0) stop("No dynamic data found.")
    yrs <- as.numeric(rownames(dyn))
    est <- dyn$Estimate
    lo  <- dyn$CI_Lower
    hi  <- dyn$CI_Upper

    if (is.null(ylim)) {
      rng  <- range(lo, hi, na.rm = TRUE)
      ylim <- c(min(rng[1], 0), max(rng[2], 0)) * 1.05
    }

    plot(yrs, est, type = "n", ylim = ylim,
         xlab = xlab %||% "Time",
         ylab = ylab %||% "Coefficients",
         main = main, mgp = c(2, .8, 0), ...)

    shade_treatment(x$tr_period, ylim)
    abline(h = 0, col = "gray50", lwd = 2, lty = 2)

    if (connected) lines(yrs, est, lwd = 2)

    ok <- lo != hi
    if (any(ok)) {
      arrows(yrs[ok], lo[ok], yrs[ok], hi[ok],
             angle = 90, code = 3, length = 0.05)
    }
    points(yrs, est, pch = 16)
    return(invisible(NULL))
  }

  ## ------------------------------------------------------------------
  ## OVERLAP
  ## ------------------------------------------------------------------
  ps <- x$ps; G <- x$G
  if (is.null(ps) || is.null(G)) stop("No propensity scores found.")

  brks <- seq(0, 1, 0.02)
  h0 <- hist(ps[G == 0], breaks = brks, plot = FALSE)
  h1 <- hist(ps[G == 1], breaks = brks, plot = FALSE)
  h0$density <- h0$counts / sum(h0$counts)
  h1$density <- h1$counts / sum(h1$counts)

  maxd <- max(h0$density, h1$density)
  ylim2 <- ylim %||% c(-maxd, maxd) * 1.1

  plot(h0, freq = FALSE, col = adjustcolor(group_colors[1], .6), ylim = ylim2,
       xlab = xlab %||% "Propensity Score",
       main = main %||% "Overlap of Propensity Scores")
  h1$density <- -h1$density
  plot(h1, freq = FALSE, add = TRUE,
       col = adjustcolor(group_colors[2], .6))
  abline(h = 0, lty = 2)
  legend("topright", legend = group_labels,
         fill = adjustcolor(group_colors, .6), bty = "n")

  invisible(NULL)
}
