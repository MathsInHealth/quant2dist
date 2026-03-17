# Greek symbol constants used in parameter display labels
.PAR_SYMS <- list(
  norm    = c("\u03bc", "\u03c3"),          # mu, sigma
  lnorm   = c("\u03bc", "\u03c3"),
  exp     = c("\u03bb"),                    # lambda
  beta    = c("\u03b1", "\u03b2"),          # alpha, beta
  gamma   = c("\u03b1", "\u03b2"),
  weibull = c("\u03bb", "\u03ba"),          # lambda, kappa
  pnorm   = c("\u03bb", "\u03bc", "\u03c3") # lambda, mu, sigma
)

.DIST_LABEL <- c(
  norm    = "normal",       lnorm   = "lognormal",
  exp     = "exponential",  beta    = "beta",
  gamma   = "gamma",        weibull = "weibull",
  pnorm   = "power normal"
)

# Format par1/par2/par3 values as "sym: val  sym: val"
.fmt_pars <- function(dist_name, row) {
  nd   <- normalize_dist(dist_name)
  syms <- .PAR_SYMS[[nd]]
  if (is.null(syms)) return("")
  pc   <- paste0("par", seq_along(syms))
  pc   <- pc[pc %in% names(row)]
  vals <- unlist(row[pc])
  ok   <- !is.na(vals)
  if (!any(ok)) return("")
  paste(paste0(syms[ok], ": ", formatC(vals[ok], format = "f", digits = 3)),
        collapse = "  ")
}

# Map any family name to a display label
.disp_name <- function(d) {
  lb <- .DIST_LABEL[normalize_dist(d)]
  if (!is.na(lb)) unname(lb) else d
}


#' Plot distribution fits interactively
#'
#' Creates an interactive \pkg{plotly} figure with three stacked panels:
#' (1) a summary of the reported input values, (2) an expanded comparison
#' table showing parameters, log-likelihood, AIC, and key quantile summaries
#' for every fitted distribution (best highlighted in blue), and (3)
#' side-by-side density and cumulative-probability curves.  Clicking a legend
#' entry toggles that distribution on both curve panels simultaneously.
#' Reported quantile values are overlaid as red dots on the CDF panel.
#'
#' Requires the \pkg{plotly} package to be installed.  If it is not available,
#' a descriptive error message is returned with installation instructions.
#'
#' @param x A \code{q2d_default} object produced by \code{\link{q2d}} with
#'   default arguments (i.e. \code{justbest = FALSE}, \code{return_df = TRUE}).
#'   The object carries the original \code{samp} input and the \code{criterion}
#'   used for fitting as attributes.
#' @param x_range Numeric vector of length 2 giving the x-axis range for
#'   the density and CDF panels.  \code{NULL} (default) uses the 0.1--99.9
#'   percentile range of the best-fitting distribution with 10 percent padding.
#' @param n_points Integer number of evaluation points for each curve.
#'   Default 400.
#'
#' @return A \code{\link[plotly]{plotly}} figure object.
#'
#' @examples
#' \dontrun{
#' # IQR + median
#' fit <- q2d(c("25\%" = 10, "50\%" = 15, "75\%" = 20, n = 100))
#' plotly_q2d(fit)
#'
#' # Five-number summary
#' fit <- q2d(c(min = 2, "25\%" = 10, "50\%" = 15, "75\%" = 20, max = 35, n = 100))
#' plotly_q2d(fit)
#'
#' # 5th, 50th, 95th percentiles, select by BIC
#' fit <- q2d(c("5\%" = 7.6, "50\%" = 7.8, "95\%" = 8.1, n = 92), criterion = "BIC")
#' plotly_q2d(fit)
#' }
#'
#' @export
plotly_q2d <- function(x, x_range = NULL, n_points = 400L) {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "The 'plotly' package is required for plotly_q2d() but is not installed.\n",
      "Install it with: install.packages(\"plotly\")",
      call. = FALSE
    )
  }

  if (!inherits(x, "q2d_default")) {
    stop(
      "'x' must be a 'q2d_default' object returned by q2d() with default ",
      "arguments (justbest = FALSE, return_df = TRUE).",
      call. = FALSE
    )
  }

  fit       <- x
  samp      <- attr(x, "samp")
  criterion <- attr(x, "criterion")
  n_points  <- as.integer(n_points)

  # ---- 1. Extract best-fit information ------------------------------------
  best_row  <- fit[fit$best == 1L, , drop = FALSE]
  best_dist <- best_row$dist[1L]
  plottable <- fit[is.finite(fit$logLik), , drop = FALSE]

  # ---- 2. x-axis range ----------------------------------------------------
  if (is.null(x_range)) {
    x_lo <- if ("0.1%"  %in% names(best_row)) as.numeric(best_row[["0.1%"]])  else as.numeric(best_row[["2.5%"]])
    x_hi <- if ("99.9%" %in% names(best_row)) as.numeric(best_row[["99.9%"]]) else as.numeric(best_row[["97.5%"]])
    if (!is.finite(x_lo) || !is.finite(x_hi) || x_lo >= x_hi) {
      x_lo <- as.numeric(best_row[["mean"]]) - 4 * as.numeric(best_row[["sd"]])
      x_hi <- as.numeric(best_row[["mean"]]) + 4 * as.numeric(best_row[["sd"]])
    }
    pad     <- (x_hi - x_lo) * 0.10
    x_range <- c(x_lo - pad, x_hi + pad)
  }
  x_seq <- seq(x_range[1L], x_range[2L], length.out = n_points)

  # ---- 3. Visual constants ------------------------------------------------
  COL_BEST  <- "#1B4F72"
  COL_OTHER <- "#AAAAAA"
  COL_OBS   <- "firebrick"

  # ---- 4. Compute density and CDF curves (needed for y-range) -------------
  n_fit     <- nrow(fit)
  is_best   <- fit$dist == best_dist
  failed    <- !is.finite(fit$logLik)
  n_plot    <- nrow(plottable)
  dens_list <- vector("list", n_plot)
  cdf_list  <- vector("list", n_plot)

  for (i in seq_len(n_plot)) {
    row_i    <- plottable[i, , drop = FALSE]
    dist_i   <- row_i$dist
    par_cols <- grep("^par[0-9]+$", names(row_i), value = TRUE)
    pars_i   <- as.numeric(unlist(row_i[par_cols]))
    pars_i   <- pars_i[!is.na(pars_i)]

    dfun <- tryCatch(get_dfun(dist_i), error = function(e) NULL)
    pfun <- tryCatch(get_pfun(dist_i), error = function(e) NULL)

    dens_y <- if (!is.null(dfun))
      tryCatch(suppressWarnings(do.call(dfun, c(list(x_seq), as.list(pars_i)))),
               error = function(e) rep(NA_real_, n_points))
    else rep(NA_real_, n_points)

    cdf_y <- if (!is.null(pfun))
      tryCatch(suppressWarnings(do.call(pfun, c(list(x_seq), as.list(pars_i)))),
               error = function(e) rep(NA_real_, n_points))
    else rep(NA_real_, n_points)

    dens_y[!is.finite(dens_y)] <- NA_real_
    cdf_y[!is.finite(cdf_y)]   <- NA_real_
    dens_list[[i]] <- dens_y
    cdf_list[[i]]  <- cdf_y
  }

  # Density y-range: start at 0, 10% headroom above max
  max_dens_val <- max(vapply(dens_list,
                             function(y) if (any(is.finite(y))) max(y, na.rm = TRUE) else 0,
                             numeric(1L)),
                      na.rm = TRUE)
  if (!is.finite(max_dens_val) || max_dens_val <= 0) max_dens_val <- 1
  dens_y_range <- c(0, max_dens_val * 1.10)

  # ---- 5. Expanded comparison table ---------------------------------------
  fmt_val <- function(col, digits = 2) {
    vals <- if (col %in% names(fit)) as.numeric(fit[[col]]) else rep(NA_real_, n_fit)
    ifelse(!is.finite(vals),
           "\u2014",
           formatC(vals, format = "f", digits = digits))
  }

  par_str <- vapply(seq_len(n_fit), function(i)
    if (failed[i]) "\u2014" else .fmt_pars(fit$dist[i], fit[i, ]),
    character(1L))

  # Distribution labels: star marks the best
  dn_str <- vapply(seq_len(n_fit), function(i) {
    lbl <- .disp_name(fit$dist[i])
    if (is_best[i]) paste0(lbl, " \u2605") else lbl
  }, character(1L))

  bg_row  <- ifelse(is_best, "#D6E4F0", "white")
  fg_row  <- ifelse(is_best, COL_BEST,  "#333333")

  comp_hdr <- list(
    values = list(
      "<b>Distribution</b>", "<b>Parameters</b>",
      "<b>logLik</b>", "<b>AIC</b>",
      "<b>Mean</b>", "<b>SD</b>",
      "<b>2.5%</b>", "<b>25%</b>", "<b>50%</b>", "<b>75%</b>", "<b>97.5%</b>"
    ),
    fill   = list(color = COL_BEST),
    font   = list(color = "white", size = 10),
    align  = c("left", "left", rep("right", 9L)),
    line   = list(color = "white", width = 1),
    height = 22
  )
  comp_cells <- list(
    values = list(
      dn_str,
      par_str,
      fmt_val("logLik",  2),
      fmt_val("AIC",     2),
      fmt_val("mean",    3),
      fmt_val("sd",      3),
      fmt_val("2.5%",    3),
      fmt_val("25%",     3),
      fmt_val("50%",     3),
      fmt_val("75%",     3),
      fmt_val("97.5%",   3)
    ),
    fill   = list(color = list(bg_row, bg_row, bg_row, bg_row, bg_row,
                               bg_row, bg_row, bg_row, bg_row, bg_row, bg_row)),
    font   = list(size = 10, color = fg_row),
    align  = c("left", "left", rep("right", 9L)),
    height = 20
  )

  # ---- 6. Input summary table ---------------------------------------------
  on_samp     <- names(samp)
  quant_names <- on_samp[on_samp != "n"]
  n_val       <- samp[["n"]]
  q_vals      <- formatC(samp[quant_names], format = "g", digits = 4)

  inp_hdr <- list(
    values = as.list(c("<b>n</b>", paste0("<b>", quant_names, "</b>"))),
    fill   = list(color = COL_BEST),
    font   = list(color = "white", size = 11),
    align  = "center",
    line   = list(color = "white", width = 1),
    height = 22
  )
  inp_cells <- list(
    values = as.list(unname(c(formatC(n_val, format = "d"), q_vals))),
    fill   = list(color = "#EBF5FB"),
    font   = list(size = 11),
    align  = "center",
    height = 20
  )

  # ---- 7. Reported quantile dots for the CDF panel -----------------------
  obs_x <- obs_p <- numeric(0)
  obs_lab <- character(0)

  pct_names <- unname(grep("%", on_samp, fixed = TRUE, value = TRUE))
  for (pn in pct_names) {
    pv <- as.numeric(sub("%", "", pn)) / 100
    if (is.finite(pv) && pv > 0 && pv < 1) {
      obs_x   <- c(obs_x, samp[[pn]])
      obs_p   <- c(obs_p, pv)
      obs_lab <- c(obs_lab, pn)
    }
  }

  # For min/max, evaluate the best-fit CDF at those x values
  pfun_b  <- tryCatch(get_pfun(best_dist), error = function(e) NULL)
  bpc     <- grep("^par[0-9]+$", names(best_row), value = TRUE)
  bpars   <- as.numeric(unlist(best_row[bpc]))
  bpars   <- bpars[!is.na(bpars)]
  eval_cdf <- function(xv) {
    if (is.null(pfun_b)) return(NA_real_)
    v <- tryCatch(suppressWarnings(
      do.call(pfun_b, c(list(xv), as.list(bpars)))), error = function(e) NA_real_)
    if (is.finite(v)) v else NA_real_
  }
  if ("min" %in% on_samp) {
    obs_x   <- c(obs_x, samp[["min"]])
    obs_p   <- c(obs_p, eval_cdf(samp[["min"]]))
    obs_lab <- c(obs_lab, "min")
  }
  if ("max" %in% on_samp) {
    obs_x   <- c(obs_x, samp[["max"]])
    obs_p   <- c(obs_p, eval_cdf(samp[["max"]]))
    obs_lab <- c(obs_lab, "max")
  }
  keep    <- !is.na(obs_p)
  obs_x   <- obs_x[keep]
  obs_p   <- obs_p[keep]
  obs_lab <- obs_lab[keep]

  # ---- 8. Domain layout (paper coords, 0 = bottom, 1 = top) --------------
  # Three stacked sections from bottom to top:
  #   [0.00, 0.43]  graphs (density + CDF side by side)
  #   [0.47, 0.80]  comparison table
  #   [0.83, 0.93]  input summary table
  # Gaps hold graph/section title annotations.
  y_graphs <- c(0.00, 0.43)
  y_ctbl   <- c(0.47, 0.80)
  y_input  <- c(0.83, 0.93)

  # x domains for graph axes (full-width, legend sits in right margin)
  x_dens   <- c(0.03, 0.46)
  x_cdf    <- c(0.52, 0.97)

  # Midpoints of the gaps — annotation anchor points
  y_gap_graphs <- (y_graphs[2L] + y_ctbl[1L]) / 2    # ~0.45
  y_gap_inp    <- (y_ctbl[2L]   + y_input[1L]) / 2   # ~0.815

  # ---- 9. Build figure ---------------------------------------------------
  fig <- plotly::plot_ly()

  # Table 1: input summary
  fig <- plotly::add_trace(fig, type = "table",
    domain  = list(x = c(0, 1), y = y_input),
    header  = inp_hdr,
    cells   = inp_cells,
    visible = TRUE)

  # Table 2: expanded comparison table
  fig <- plotly::add_trace(fig, type = "table",
    domain      = list(x = c(0, 1), y = y_ctbl),
    header      = comp_hdr,
    cells       = comp_cells,
    columnwidth = list(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    visible     = TRUE)

  # Distribution curves (density + CDF paired per legend group)
  for (i in seq_len(n_plot)) {
    row_i     <- plottable[i, , drop = FALSE]
    dist_i    <- row_i$dist
    disp_i    <- .disp_name(dist_i)
    is_best_i <- dist_i == best_dist
    clr_i     <- if (is_best_i) COL_BEST else COL_OTHER
    lwd_i     <- if (is_best_i) 2.5 else 1.2

    dens_y <- dens_list[[i]]
    cdf_y  <- cdf_list[[i]]

    # Label position for density: peak of the curve
    peak_idx <- if (any(is.finite(dens_y))) which.max(dens_y) else ceiling(n_points / 2L)
    lab_dx   <- x_seq[peak_idx]
    lab_dy   <- dens_y[peak_idx]

    # Label position for CDF: stagger across distributions to reduce overlap
    cdf_target <- 0.50 + 0.07 * ((i - 1L) %% 5L)
    cdf_idx    <- which.min(abs(cdf_y - cdf_target))
    lab_cx     <- x_seq[cdf_idx]
    lab_cy     <- cdf_y[cdf_idx]

    # Density line
    fig <- plotly::add_trace(fig,
      x = x_seq, y = dens_y,
      type        = "scatter", mode = "lines",
      xaxis       = "x",       yaxis = "y",
      name        = disp_i,
      legendgroup = disp_i,
      showlegend  = TRUE,
      line        = list(color = clr_i, width = lwd_i),
      hovertemplate = paste0("<b>", disp_i,
        "</b><br>x: %{x:.3f}<br>density: %{y:.4f}<extra></extra>"))

    # Density label
    if (is.finite(lab_dy)) {
      fig <- plotly::add_trace(fig,
        x = lab_dx, y = lab_dy,
        type         = "scatter", mode = "text",
        xaxis        = "x",       yaxis = "y",
        text         = disp_i,
        textposition = "top center",
        textfont     = list(color = clr_i,
                            size  = if (is_best_i) 11L else 9L),
        legendgroup  = disp_i,
        showlegend   = FALSE,
        hoverinfo    = "skip")
    }

    # CDF line
    fig <- plotly::add_trace(fig,
      x = x_seq, y = cdf_y,
      type        = "scatter", mode = "lines",
      xaxis       = "x2",      yaxis = "y2",
      name        = disp_i,
      legendgroup = disp_i,
      showlegend  = FALSE,
      line        = list(color = clr_i, width = lwd_i),
      hovertemplate = paste0("<b>", disp_i,
        "</b><br>x: %{x:.3f}<br>P: %{y:.3f}<extra></extra>"))

    # CDF label
    if (is.finite(lab_cy)) {
      fig <- plotly::add_trace(fig,
        x = lab_cx, y = lab_cy,
        type         = "scatter", mode = "text",
        xaxis        = "x2",      yaxis = "y2",
        text         = disp_i,
        textposition = "top right",
        textfont     = list(color = clr_i,
                            size  = if (is_best_i) 11L else 9L),
        legendgroup  = disp_i,
        showlegend   = FALSE,
        hoverinfo    = "skip")
    }
  }

  # Reported quantile dots on CDF
  if (length(obs_x) > 0L) {
    fig <- plotly::add_trace(fig,
      x    = obs_x,  y = obs_p,
      type = "scatter", mode = "markers+text",
      xaxis = "x2",  yaxis = "y2",
      name  = "Reported",
      marker = list(color = COL_OBS, size = 8L, symbol = "circle"),
      text  = obs_lab,
      textposition = "top right",
      textfont = list(color = COL_OBS, size = 9L),
      showlegend   = TRUE,
      hovertemplate = paste0(
        "<b>%{text}</b><br>x: %{x:.3f}<br>P: %{y:.3f}<extra></extra>"))
  }

  # ---- 10. Layout --------------------------------------------------------
  fig <- plotly::layout(fig,
    paper_bgcolor = "white",
    plot_bgcolor  = "white",
    margin = list(l = 10, r = 130, t = 30, b = 10),
    xaxis = list(
      domain    = x_dens,
      anchor    = "y",
      title     = "",
      showgrid  = TRUE,
      gridcolor = "#EEEEEE",
      range     = x_range
    ),
    yaxis = list(
      domain    = y_graphs,
      anchor    = "x",
      title     = "Density",
      showgrid  = TRUE,
      gridcolor = "#EEEEEE",
      range     = dens_y_range   # explicit [0, max*1.1] — aligns zero with CDF
    ),
    xaxis2 = list(
      domain    = x_cdf,
      anchor    = "y2",
      title     = "",
      showgrid  = TRUE,
      gridcolor = "#EEEEEE",
      range     = x_range
    ),
    yaxis2 = list(
      domain    = y_graphs,          # same domain as yaxis — zeros align
      anchor    = "x2",
      title     = "Cumulative probability",
      showgrid  = TRUE,
      gridcolor = "#EEEEEE",
      range     = c(0, 1)
    ),
    legend = list(
      orientation = "v",
      x           = 1.01,
      y           = 0.22,
      xanchor     = "left",
      bgcolor     = "white",
      bordercolor = "#CCCCCC",
      borderwidth = 1
    ),
    annotations = list(
      # "Reported quantile information" — above input table
      list(
        x = 0.5, y = (y_input[2L] + 1.00) / 2,
        xref = "paper", yref = "paper",
        text = "<b>Reported quantile information</b>",
        showarrow = FALSE,
        font = list(size = 12, color = COL_BEST),
        xanchor = "center", yanchor = "middle"
      ),
      # "Density" — in the gap between comparison table and graph panels
      list(
        x = mean(x_dens), y = y_gap_graphs,
        xref = "paper", yref = "paper",
        text = "<b>Density</b>",
        showarrow = FALSE,
        font = list(size = 13, color = COL_BEST),
        xanchor = "center", yanchor = "middle"
      ),
      # "Cumulative probability" — in the same gap, right panel
      list(
        x = mean(x_cdf), y = y_gap_graphs,
        xref = "paper", yref = "paper",
        text = "<b>Cumulative probability</b>",
        showarrow = FALSE,
        font = list(size = 13, color = COL_BEST),
        xanchor = "center", yanchor = "middle"
      )
    )
  )

  fig
}
