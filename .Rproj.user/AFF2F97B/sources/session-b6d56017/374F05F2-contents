#' Plot Results from FDID Analysis
#'
#' Provides visualizations for FDID results, including raw means, dynamic effects, propensity score overlap, and method comparisons.
#'
#' @param x An `fdid` object or a list of multiple `fdid` objects.
#' @param type Type of plot to produce: `"raw"`, `"dynamic"`, `"overlap"`, or `"est"`.
#' @param alpha_shade Transparency for shading the treatment period.
#' @param col_treated,col_control Colors for treated and control groups.
#' @param xlab,ylab,main Axis labels and main title.
#' @param ylim Y-axis limits. Default is `NULL`.
#'
#' @return Produces a plot. Nothing is returned.
#'
#' @examples
#' plot.fdid(fdid_results, type = "raw")
#' plot.fdid(fdid_results, type = "dynamic")
#'
#' @export
plot.fdid <- function(x,
                      type = c("raw", "dynamic", "overlap", "est"),
                      alpha_shade = 0.2,
                      col_treated = "red",
                      col_control = "blue",
                      xlab = NULL,
                      ylab = NULL,
                      main = NULL,
                      ylim = NULL) {
  type <- match.arg(type)

  # -------------------------------------------------------------------
  # Helper for shading the treatment period, if it exists
  # -------------------------------------------------------------------
  shade_treatment <- function(tr_period, yrange, col = "gray") {
    tmin <- min(tr_period)
    tmax <- max(tr_period)
    rect(xleft = tmin-0.5,
         ybottom = yrange[1],
         xright = tmax+0.5,
         ytop = yrange[2],
         col = adjustcolor(col, alpha.f = alpha_shade),
         border = NA)
  }

  # -------------------------------------------------------------------
  # 1) Plot RAW means
  # -------------------------------------------------------------------
  if (type == "raw") {
    if (is.null(x$raw_means)) {
      stop("No 'raw_means' found in 'fdid' object.")
    }
    rawdf <- x$raw_means  # data.frame(year, group, meanY)
    # Quick checks
    if (!all(c("year","group","meanY") %in% names(rawdf))) {
      stop("raw_means missing required columns (year, group, meanY).")
    }
    all_years <- sort(unique(rawdf$year))
    if (is.null(ylim)) {
      ylim <- range(rawdf$meanY, na.rm=TRUE)
    }
    if (is.null(xlab)) xlab <- "Year"
    if (is.null(ylab)) ylab <- "Mean outcome"

    plot(NULL, xlim = range(all_years), ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, ...)
    shade_treatment(x$tr_period, ylim)

    # Plot lines: Treated vs Control
    treated_sub  <- rawdf[rawdf$group=="Treated", ]
    control_sub  <- rawdf[rawdf$group=="Control", ]

    lines(treated_sub$year, treated_sub$meanY, col=col_treated, lwd=2)
    points(treated_sub$year, treated_sub$meanY, pch=16, col=col_treated)

    lines(control_sub$year, control_sub$meanY, col=col_control, lty=2, lwd=2)
    points(control_sub$year, control_sub$meanY, pch=1, col=col_control)

    legend("topleft",
           legend = c("Treated","Control"),
           col = c(col_treated,col_control),
           pch = c(16,1), lty = c(1,2), bty="n")

    return(invisible(NULL))
  }

  # -------------------------------------------------------------------
  # 2) Plot DYNAMIC estimates (year-by-year FDID)
  # -------------------------------------------------------------------
  if (type == "dynamic") {
    if (is.null(x$dynamic) || nrow(x$dynamic)==0) {
      stop("No dynamic data found in 'x$dynamic'.")
    }
    dyn <- x$dynamic
    yrs <- as.numeric(rownames(dyn))
    est <- dyn$Estimate
    lo  <- dyn$CI_Lower
    hi  <- dyn$CI_Upper

    if (is.null(ylim)) {
      ylim <- range(lo, hi, na.rm=TRUE)
    }
    if (is.null(xlab)) xlab <- "Year"
    if (is.null(ylab)) ylab <- "Treatment Effect (Dynamic FDID)"

    plot(yrs, est, type="n", ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...)
    shade_treatment(x$tr_period, ylim)
    abline(h=0, col="gray50", lwd=2, lty=2)

    # error bars
    arrows(yrs, lo, yrs, hi, angle=90, code=3, length=0.05, col="black")
    points(yrs, est, pch=16)

    return(invisible(NULL))
  }

  # -------------------------------------------------------------------
  # 3) Plot OVERLAP of propensity scores
  # -------------------------------------------------------------------
  if (type == "overlap") {
    if (is.null(x$ps) || is.null(x$G)) {
      stop("'fdid' object has no stored propensity scores (x$ps) or no x$G.")
    }
    ps <- x$ps
    G  <- x$G
    # Very simple mirrored histogram approach:

    brks <- seq(0, 1, by=0.02)
    # control group histogram
    hCo <- hist(ps[G==0], breaks=brks, plot=FALSE)
    # treated group histogram
    hTr <- hist(ps[G==1], breaks=brks, plot=FALSE)

    # convert counts to fraction (density-like)
    hCo$density <- hCo$counts / sum(hCo$counts, na.rm=TRUE)
    hTr$density <- hTr$counts / sum(hTr$counts, na.rm=TRUE)
    ymax <- max(hCo$density, hTr$density)

    if (is.null(ylim)) {
      ylim <- c(-ymax, ymax)
    }
    if (is.null(xlab)) xlab <- "Propensity Score"
    if (is.null(main)) main <- "Overlap of Propensity Scores"

    plot(hCo, freq=FALSE, col=adjustcolor("blue", 0.4),
         ylim=ylim, main=main, xlab=xlab, border=NA)
    # Flip the treated histogram to negative for a mirrored effect
    hTr$density <- -hTr$density
    plot(hTr, freq=FALSE, col=adjustcolor("red", 0.4),
         add=TRUE, border=NA)
    abline(h=0, lty=2, col="gray50")

    legend("topright", legend=c("Control","Treated"),
           fill=c("blue","red"), border=NA, bty="n", cex=0.9)
    return(invisible(NULL))
  }

  # -------------------------------------------------------------------
  # 4) Plot COMPARISON of ESTIMATES if multiple methods are used
  # -------------------------------------------------------------------
  if (type == "est") {
    # We assume 'x' is actually a *list* of fdid objects
    if (!is.list(x)) {
      stop("For type='est', supply a list of 'fdid' objects, each from a different method.")
    }
    # each element should have x[[i]]$static => single row
    est_list <- lapply(x, function(obj) {
      st <- obj$static
      if (is.null(st) || nrow(st)==0) {
        stop("One of the list elements has no $static result.")
      }
      data.frame(method   = obj$method,
                 Estimate = st$Estimate,
                 CI_Lower = st$CI_Lower,
                 CI_Upper = st$CI_Upper)
    })
    est_df <- do.call(rbind, est_list)
    n <- nrow(est_df)
    if (is.null(ylim)) {
      ylim <- range(c(est_df$CI_Lower, est_df$CI_Upper), na.rm=TRUE)
    }
    if (is.null(xlab)) xlab <- "Estimate"
    if (is.null(ylab)) ylab <- ""
    if (is.null(main)) main <- "Comparison of FDID Estimates"

    # base R dot-whisker style
    par(mar = c(5, 8, 2, 2))
    plot(NA, xlim = ylim, ylim = c(0.5, n + 0.5),
         ylab = "", xlab = xlab, main = main, yaxt="n")
    abline(v = 0, lty=2, col="gray50") # reference line
    axis(2, at = seq_len(n), labels = est_df$method, las=1)

    for (i in seq_len(n)) {
      segments(x0 = est_df$CI_Lower[i], x1 = est_df$CI_Upper[i],
               y0 = i, y1 = i, lwd=2)
    }
    points(est_df$Estimate, seq_len(n), pch=16)
    box()

    return(invisible(NULL))
  }
}
