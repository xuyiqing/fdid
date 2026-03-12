#' Plot Results from FDID Analysis
#'
#' Provides visualisations for FDID results, including raw means, dynamic
#' effects, and propensity-score overlap.  The comparison plot of multiple
#' methods has been removed; use \code{plot.fdid_list()} for that.
#'
#' @param x   An \code{fdid} object.
#' @param type One of \code{"raw"}, \code{"dynamic"}, or \code{"overlap"}.
#' @param connected Logical; if \code{TRUE}, connects points with lines in
#'        the "raw" and "dynamic" plots.  Default is \code{FALSE}.
#' @param ci Logical; if \code{TRUE}, draw 95\% CIs when available. Default is \code{TRUE}.
#' @param shade_periods Shaded intervals on the time axis. Default uses \code{x$tr_period}, i.e. event periods.
#'        Set to \code{NULL} to remove shaded area.
#' @param alpha_shade Transparency for shading the treatment period.
#' @param palette A palette name from \strong{RColorBrewer}.  Default \code{"Set2"}.
#' @param group_labels Labels for the two groups.
#' @param xlab,ylab,main Axis labels and main title.
#' @param ylim Y-axis limits.  Default \code{NULL} (computed automatically).
#' @param ...  Additional graphics parameters.
#'
#' @return Produces a plot; invisibly returns \code{NULL}.
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
#' plot(result, type = "raw")
#' plot(result, type = "dynamic")
#' }
#' @author Rivka Lipkovitz, Enhan Liu
#' @importFrom rlang %||%
#' @export
plot.fdid <- function(x,
                      type = c("raw", "dynamic", "overlap"),
                      connected = FALSE,
                      ci = TRUE,
                      shade_periods = x$tr_period,
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
    
    # remove shading: user passes NULL/empty
    if (is.null(tr_period) || length(tr_period) == 0 ||
        !is.finite(alpha_shade) || alpha_shade <= 0) {
      return(invisible(NULL))
    }
    
    draw_one <- function(tp) {
      tp <- as.numeric(tp)
      tp <- tp[!is.na(tp)]
      if (length(tp) == 0) return(invisible(NULL))
      
      rect(min(tp) - 0.5,
           yrange[1] - 10,
           max(tp) + 0.5,
           yrange[2] + 10,
           col = adjustcolor(col, alpha.f = alpha_shade),
           border = NA)
      invisible(NULL)
    }
    
    if (is.list(tr_period)) {
      for (tp in tr_period) draw_one(tp)
    } else {
      draw_one(tr_period)
    }
    
    invisible(NULL)
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
    
    # If CI columns exist (or seY exists), prepare 95% CI for plotting
    has_ci <- ci && all(c("CI_Lower", "CI_Upper") %in% names(rawdf))
    
    if (is.null(ylim)) {
      if (has_ci) {
        tmp <- range(rawdf$CI_Lower, rawdf$CI_Upper, na.rm = TRUE)
      } else {
        tmp <- range(rawdf$meanY, na.rm = TRUE)
      }
      ylim <- c(min(tmp[1], 0.99 * tmp[1]), max(tmp[2], 1.01 * tmp[2]))
    }

    plot(NULL, xlim = range(times), ylim = ylim,
         xlab = xlab %||% "Time",
         ylab = ylab %||% "Mean outcome",
         main = main, mgp = c(2, .8, 0), ...)
    shade_treatment(shade_periods, ylim)

    g1 <- rawdf[rawdf$group == "Group 1", ]
    g0 <- rawdf[rawdf$group == "Group 0", ]

    if (connected) {
      lines(g1$time, g1$meanY, col = group_colors[2], lwd = 2)
      lines(g0$time, g0$meanY, col = group_colors[1], lwd = 2)
    }
    
    # 95% CI error bars (default on when available)
    if (exists("has_ci") && isTRUE(has_ci)) {
      ok1 <- !is.na(g1$CI_Lower) & !is.na(g1$CI_Upper)
      ok0 <- !is.na(g0$CI_Lower) & !is.na(g0$CI_Upper)
      if (any(ok1)) {
        arrows(g1$time[ok1], g1$CI_Lower[ok1], g1$time[ok1], g1$CI_Upper[ok1],
               angle = 90, code = 3, length = 0.05, col = group_colors[2])
      }
      if (any(ok0)) {
        arrows(g0$time[ok0], g0$CI_Lower[ok0], g0$time[ok0], g0$CI_Upper[ok0],
               angle = 90, code = 3, length = 0.05, col = group_colors[1])
      }
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
      if (isTRUE(ci)) {
        rng  <- range(lo, hi, na.rm = TRUE)
      } else {
        rng  <- range(est, na.rm = TRUE)
      }
      ylim <- c(min(rng[1], 0), max(rng[2], 0)) * 1.05
    }

    plot(yrs, est, type = "n", ylim = ylim,
         xlab = xlab %||% "Time",
         ylab = ylab %||% "Coefficients",
         main = main, mgp = c(2, .8, 0), ...)

    shade_treatment(shade_periods, ylim)
    abline(h = 0, col = "gray50", lwd = 2, lty = 2)

    if (connected) lines(yrs, est, lwd = 2)

    if (isTRUE(ci)) {
      ok <- lo != hi
      if (any(ok)) {
        arrows(yrs[ok], lo[ok], yrs[ok], hi[ok],
               angle = 90, code = 3, length = 0.05)
      }
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
