#' Compute fallback starting values from quantiles
#'
#' Uses method-of-moments approximations (Luo et al. 2018, Wan et al. 2014)
#' to derive starting parameter values when \code{estmeansd::qe.fit} fails.
#'
#' @param quants Numeric vector of quantile values.
#' @param n Sample size.
#' @param scenario One of \code{"S1"}, \code{"S2"}, \code{"S3"}.
#' @return A list of named starting values for each distribution family.
#'
#' @importFrom stats qnorm
#' @keywords internal
fallback_start_values <- function(quants, n, scenario) {
  if(scenario == "S3") {
    # 5-number summary: min, q1, median, q3, max
    a <- quants[1]; q1 <- quants[2]; m <- quants[3]; q3 <- quants[4]; b <- quants[5]
    mean_hat <- (a + 2*q1 + 2*m + 2*q3 + b) / 8
    sd_hat <- (b - a) / (4 * qnorm((n - 0.375) / (n + 0.25)))
  } else if(scenario == "S1") {
    # 3-number summary: min, median, max
    a <- quants[1]; m <- quants[2]; b <- quants[3]
    mean_hat <- (a + 2*m + b) / 4
    sd_hat <- (b - a) / (2 * qnorm((n - 0.375) / (n + 0.25)))
  } else {
    # S2: q1, median, q3
    q1 <- quants[1]; m <- quants[2]; q3 <- quants[3]
    mean_hat <- (q1 + m + q3) / 3
    sd_hat <- (q3 - q1) / (2 * qnorm((0.75 * n - 0.125) / (n + 0.25)))
  }

  k_start <- (mean_hat / sd_hat)^1.086
  list(
    norm.mu.start    = mean_hat,
    norm.sigma.start = sd_hat,
    lnorm.mu.start    = log(mean_hat / sqrt(1 + (sd_hat / mean_hat)^2)),
    lnorm.sigma.start = sqrt(log(1 + (sd_hat / mean_hat)^2)),
    gamma.shape.start = mean_hat^2 / sd_hat^2,
    gamma.rate.start  = mean_hat / sd_hat^2,
    beta.shape1.start = mean_hat * (((mean_hat * (1 - mean_hat)) / (sd_hat^2)) - 1),
    beta.shape2.start = (1 - mean_hat) * (((mean_hat * (1 - mean_hat)) / (sd_hat^2)) - 1),
    weibull.shape.start = k_start,
    weibull.scale.start = mean_hat / gamma((k_start + 1) / k_start),
    pnorm.lambda.start = 1,
    pnorm.mu.start    = mean_hat,
    pnorm.sigma.start = sd_hat
  )
}

#' Estimate starting values for distribution families
#'
#' Prepares starting parameter values for optimization by combining
#' \code{estmeansd::qe.fit} results with method-of-moments fallbacks.
#'
#' @param samp Named numeric vector with quantile values and \code{n}.
#' @param families Character vector of distribution family names to return
#'   starting values for. Defaults to the six built-in families. Unknown
#'   families fall back to normal-distribution starting values.
#' @param retorig Logical; if \code{TRUE}, return raw \code{qe.fit} output.
#' @param print_debug Logical; if \code{TRUE}, print debug information.
#' @return A named list of starting parameter vectors, one per distribution.
#' @keywords internal
estimate_start_values <- function(samp,
                                  families = c("normal", "lognormal", "exponential",
                                               "beta", "gamma", "weibull"),
                                  retorig = FALSE, print_debug = FALSE) {
  stopifnot('n' %in% names(samp))
  thisn <- samp[['n']]

  on <- names(samp)
  if("sd" %in% on && "mean" %in% on) samp <- c("25%" = samp[['mean']] - 0.675* samp[['sd']], '50%' = samp[['mean']], "75%" = samp[['mean']] + 0.675* samp[['sd']], 'n' = samp[['n']])
  if("se" %in% on && "mean" %in% on) samp <- c("25%" = samp[['mean']] - 0.675* samp[['se']] * sqrt(samp[['n']]), '50%' = samp[['mean']], "75%" = samp[['mean']] + 0.675* samp[['se']] * sqrt(samp[['n']]), 'n' = samp[['n']])
  stopifnot(length(these_quants <- c(grep('%', names(samp), fixed = TRUE), which(names(samp) %in% c('min', 'mean', 'max'))))>2)
  orig_samp <- samp
  on <- names(samp)
  samp <- samp[which(names(samp) %in% c("min", "25%", "50%", "75%", "max", "n"))]
  if(!"50%" %in% names(samp)) {
    if("mean" %in% on) samp <- c(samp, '50%' = orig_samp[['mean']])
  }
  if(length(samp) < 4) {
    samp <- orig_samp[these_quants]

    if(!"50%" %in% names(samp)) {
      if("mean" %in% on) samp <- {
        c(samp, '50%' = orig_samp[['mean']])
        samp <- samp[-which(names(samp) == 'mean')]
      }
    }
    if('min' %in% names(samp)) names(samp)[which(names(samp)=='min')] <- '0%'
    if('max' %in% names(samp)) names(samp)[which(names(samp)=='max')] <- '100%'


    samp_n <- samp[grep('%', names(samp), fixed = TRUE)]

    stopifnot(length(samp_n)>1)

    orig_v <- sapply(regmatches(names(samp_n), gregexpr("^\\d*\\.?\\d+(?=%)",  names(samp_n), perl = TRUE)), as.numeric)/100
    samp_n <- samp_n[order(orig_v)] # Sort
    orig_v <- sort(orig_v) # Sort


    orig_intervals <- findInterval(orig_v, (c(0, 1, 3, 5, 7, 8))/8)
    names(samp_n) <- switch(paste(orig_intervals, collapse = ""),
                            "135" = c('min', '50%', 'max'),
                            "12345" = c('min', '25%', '50%', '75%', 'max'),
                            "234" = c('25%', '50%', '75%'),
                            "134" = c('25%', '50%', '75%'),
                            "235" = c('25%', '50%', '75%'),
                            {
                              samp_n <- samp_n[c(1, ceiling(length(samp_n)/2), length(samp_n))]
                              c('25%', '50%', '75%')
                            }
    )


    samp <- c(samp_n, 'n' = thisn)

  }

  names(samp) <- c("min.val", "q1.val", "med.val", "q3.val", "max.val", "n")[match(names(samp), c("min", "25%", "50%", "75%", "max", "n"))]
  scenario <- ifelse(length(samp) == 6, "S3", ifelse("min.val" %in% names(samp), "S1", "S2"))
  quants <- samp[-which(names(samp) == "n")]
  samp <- as.list(samp)
  sampn <- samp
  sampn[['n']] <- max(3, sampn[['n']])
  if(print_debug) {
    print('quants:\r')
    print(quants)
    print('\rsampn:\r')
    print(sampn)
  }
  tmps <- fallback_start_values(quants = quants, n = sampn[['n']], scenario = scenario)
  tmp <- as.list(do.call(estmeansd::qe.fit, sampn))
  if(retorig) return(tmp)

  norm_start <- if(is.na(tmp$norm.par[1])) tmps$norm.mu.start[[1]] else tmp$norm.par[['mu']]
  sd_start   <- if(is.na(tmp$norm.par[1])) tmps$norm.sigma.start[[1]] else tmp$norm.par[['sigma']]

  get_base_stvs <- function(nf) {
    switch(nf,
      norm = c(mean = norm_start, sd = sd_start),
      lnorm = if(is.na(tmp$lnorm.par[1])) {
        c(meanlog = tmps$lnorm.mu.start[[1]], sdlog = tmps$lnorm.sigma.start[[1]])
      } else {
        c(meanlog = tmp$lnorm.par[['mu']], sdlog = tmp$lnorm.par[['sigma']])
      },
      exp = c(rate = 0.5),
      beta = if(is.na(tmp$beta.par[1])) {
        c(shape1 = tmps$beta.shape1.start[[1]], shape2 = tmps$beta.shape2.start[[1]])
      } else {
        c(tmp$beta.par)
      },
      scalebeta = if(is.na(tmp$beta.par[1])) {
        c(shape1 = tmps$beta.shape1.start[[1]], shape2 = tmps$beta.shape2.start[[1]], SCALE = 1)
      } else {
        c(tmp$beta.par, SCALE = 1)
      },
      gamma = if(is.na(tmp$gamma.par[1])) {
        c(rate = tmps$gamma.rate.start[[1]], shape = tmps$gamma.shape.start[[1]])
      } else {
        tmp$gamma.par
      },
      weibull = if(is.na(tmp$weibull.par[1])) {
        c(shape = tmps$weibull.shape.start[[1]], scale = tmps$weibull.scale.start[[1]])
      } else {
        tmp$weibull.par
      },
      # Power Normal (powdist): lambda first, then mu, sigma (matching qpnorm argument order)
      pnorm = c(lambda = tmps$pnorm.lambda.start[[1]],
                mu     = tmps$pnorm.mu.start[[1]],
                sigma  = tmps$pnorm.sigma.start[[1]]),
      # Generic fallback for unknown distributions: use normal starting values
      c(par1 = norm_start, par2 = sd_start)
    )
  }

  setNames(lapply(families, function(f) get_base_stvs(normalize_dist(f))), families)
}


#' Fit parametric distributions to reported quantile data
#'
#' Fits one or more parametric distribution families to quantile-reported
#' summary statistics using maximum likelihood estimation, and selects the
#' best-fitting distribution.
#'
#' @param samp A named numeric vector of quantile values. Names should be
#'   quantile labels such as \code{"25\%"}, \code{"50\%"}, \code{"min"},
#'   \code{"max"}, and must include \code{"n"} for sample size. May also
#'   include \code{"mean"}, \code{"sd"}, or \code{"se"}.
#' @param families Character vector of distribution family names to fit.
#'   Defaults to \code{NULL}, which fits the seven built-in families: normal,
#'   lognormal, exponential, beta, gamma, Weibull, and Power Normal
#'   (\code{"pnorm"} from the \pkg{powdist} package).
#'   Pass a character vector to restrict the set, e.g.
#'   \code{families = c("normal", "lognormal")} or
#'   \code{families = c("normal", "pnorm")}.
#'   When any quantile value is negative, families that are defined only on
#'   non-negative values (lognormal, exponential, beta, gamma, Weibull) are
#'   automatically excluded and a warning is issued; only \code{"normal"} and
#'   \code{"pnorm"} are retained.  Explicitly passing \code{families} does not
#'   suppress this filter, but you can supply only those two families directly
#'   to avoid the warning.
#' @param justbest Logical; if \code{FALSE} (default), return results for all
#'   fitted distributions. If \code{TRUE}, return only the best-fitting
#'   distribution.
#' @param hessian Logical; if \code{TRUE}, compute the Hessian matrix at the
#'   optimum for each distribution.
#' @param tol Integer; number of decimal places for comparing log-likelihoods
#'   when \code{criterion = "logLik"} (default 8).
#' @param return_df Logical; if \code{TRUE} (default), return results as a
#'   data frame.
#' @param criterion Character string specifying the information criterion used
#'   to select the best-fitting distribution. One of \code{"AIC"} (default),
#'   \code{"BIC"}, or \code{"logLik"}. AIC and BIC are minimised; logLik is
#'   maximised.
#'
#' @return By default (\code{justbest = FALSE}, \code{return_df = TRUE}), a
#'   data frame with one row per distribution, including \code{logLik},
#'   \code{AIC}, \code{BIC}, and \code{best} columns (where \code{best} ranks
#'   distributions from 1 = best to n = worst according to the chosen
#'   criterion). This object has class \code{c("q2d_default", "data.frame")}
#'   and carries \code{samp} and \code{criterion} attributes storing the
#'   original input and the criterion used. Passing \code{justbest = TRUE}
#'   returns only the best row; \code{return_df = FALSE} returns a list.
#'
#' @examples
#' # Fit to IQR + median (default: all distributions, data frame output)
#' q2d(c("25%" = 10, "50%" = 15, "75%" = 20, n = 100))
#'
#' # Fit to 5-number summary
#' q2d(c(min = 2, "25%" = 10, "50%" = 15, "75%" = 20, max = 35, n = 100))
#'
#' # Return only best-fitting distribution
#' q2d(c("25%" = 10, "50%" = 15, "75%" = 20, n = 100), justbest = TRUE)
#'
#' # Fit only normal and lognormal
#' q2d(c("25%" = 10, "50%" = 15, "75%" = 20, n = 100),
#'     families = c("normal", "lognormal"))
#'
#' @export
q2d <- function(samp, families = NULL, justbest = FALSE, hessian = FALSE,
                tol = 8, return_df = TRUE, criterion = "AIC") {
  criterion <- match.arg(criterion, choices = c("AIC", "BIC", "logLik"))

  if (is.null(families)) {
    families <- c("normal", "lognormal", "exponential", "beta", "gamma", "weibull", "pnorm")
  }

  # When no quantile data are present (only mean/sd/se), AIC and BIC penalise
  # distributions with more parameters even though the extra parameters are
  # unconstrained by the data.  This causes simpler families (e.g. exponential,
  # 1 parameter) to be selected over normal (2 parameters) regardless of fit.
  # In this case, force selection by log-likelihood and warn the user.
  non_quantile_names <- c("n", "mean", "sd", "se")
  has_quantiles <- any(!names(samp) %in% non_quantile_names)
  if (!has_quantiles && criterion %in% c("AIC", "BIC")) {
    warning(
      "No quantile data found in 'samp' (only mean/sd/se provided). ",
      "AIC/BIC model selection is unreliable without quantile constraints ",
      "because parameter count differences dominate the criterion. ",
      "Switching to criterion = \"logLik\" for best-distribution selection."
    )
    criterion <- "logLik"
  }

  # Distributions defined only on non-negative (or (0,1)) values cannot be
  # fitted when the data contain negative quantiles.  Retain only families
  # whose support covers the full real line: normal and Power Normal (pnorm).
  quant_vals <- samp[!names(samp) %in% c("n", "mean", "sd", "se")]
  if (any(quant_vals < 0, na.rm = TRUE)) {
    real_line_canonical <- c("norm", "pnorm")
    keep <- families[vapply(families, normalize_dist, character(1L)) %in%
                       real_line_canonical]
    if (length(keep) == 0L) keep <- c("normal", "pnorm")
    if (!setequal(keep, families)) {
      warning(
        "Negative quantile values detected. Only fitting distributions that ",
        "support the full real line (normal, Power Normal). ",
        "Distributions requiring non-negative values (lognormal, exponential, ",
        "beta, gamma, Weibull) have been excluded. ",
        "Use the 'families' argument to override."
      )
      families <- keep
    }
  }

  stvs <- estimate_start_values(samp, families = families)

  tmpout <- suppressWarnings(setNames(
    lapply(families, function(f) {
      fit_dist(vals = samp, par = stvs[[f]], dist = normalize_dist(f), hessian = hessian)
    }),
    families
  ))

  tn <- families
  names(tn) <- tn

  # Determine max number of parameters across all fitted families (for consistent output)
  max_pars <- max(sapply(tmpout, function(f) length(f$par)))

  tout <- lapply(tn, function(x) {
    fit_pars <- tmpout[[x]]$par
    n_pars   <- length(fit_pars)
    par_vec  <- setNames(unname(fit_pars), paste0("par", seq_along(fit_pars)))
    # Pad with NA to ensure consistent column count across families
    if (n_pars < max_pars) {
      par_vec <- c(par_vec,
                   setNames(rep(NA_real_, max_pars - n_pars),
                            paste0("par", (n_pars + 1L):max_pars)))
    }
    # +Inf logLik signals a convergence failure (optimizer returned -Inf nll);
    # treat it as -Inf so it is never selected as the best distribution.
    loglik_val <- -tmpout[[x]]$value
    if (is.infinite(loglik_val) && loglik_val > 0) loglik_val <- -Inf
    suppressWarnings(c(
      dist_moments(dist = normalize_dist(x), pars = fit_pars,
                   quantiles = c("0.1%" = 0.001, "2.5%" = 0.025, "5%" = 0.05,
                                 "25%" = 0.25, "50%" = 0.5, "75%" = 0.75,
                                 "95%" = 0.95, "97.5%" = 0.975, "99.9%" = 0.999)),
      value = loglik_val,
      par_vec
    ))
  })

  values <- sapply(tout, function(x) x[["value"]])
  n_params <- sapply(families, function(f) length(tmpout[[f]]$par))
  n_obs    <- samp[['n']]
  aics <- setNames(2 * n_params - 2 * values, families)
  bics <- setNames(n_params * log(n_obs) - 2 * values, families)

  # Add AIC and BIC to each distribution's result vector
  for (f in families) {
    tout[[f]] <- c(tout[[f]], AIC = aics[[f]], BIC = bics[[f]])
  }

  # Select best distribution based on criterion
  minv <- switch(criterion,
    logLik = {
      minvs <- round(values - max(values, na.rm = TRUE), tol)
      idx <- match(0, minvs)
      setNames(idx, names(minvs)[idx])
    },
    AIC = which.min(aics),
    BIC = which.min(bics)
  )

  tout$best <- c(tout[[minv]], "dist" = minv)
  tout$pars <- lapply(tn, function(x) c(tmpout[[x]]$par))
  tout$values <- values
  tout$n <- samp[['n']]

  if (return_df) {
    tmp <- suppressWarnings(do.call(rbind, lapply(tout[families], function(v) {
      as.data.frame(as.list(v), stringsAsFactors = FALSE, check.names = FALSE)
    })))
    tmp$logLik <- tmp$value
    tmp$dist <- rownames(tmp)
    # Rank distributions: 1 = best, n = worst, based on chosen criterion
    rank_vals <- switch(criterion,
      AIC    = rank(aics,    ties.method = "first"),
      BIC    = rank(bics,    ties.method = "first"),
      logLik = {
        # Use the same tolerance rounding as the minv selector: differences
        # smaller than 10^{-tol} from the maximum are treated as tied.
        # This prevents floating-point noise from breaking ties when all
        # distributions fit equally well (e.g. only mean/SD provided).
        rel <- round(values - max(values, na.rm = TRUE), tol)
        rank(-rel, ties.method = "first")
      }
    )
    tmp$best <- as.integer(rank_vals)
    par_cols <- sort(grep("^par[0-9]+$", names(tmp), value = TRUE))
    fixed_cols <- c("dist", par_cols, "logLik", "AIC", "BIC", "best",
                    "mean", "var", "sd", "median", "mode", "max",
                    "0.1%", "2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%", "99.9%")
    tout <- tmp[, fixed_cols[fixed_cols %in% names(tmp)]]
    if (justbest) tout <- tout[minv, ]
    if (!justbest) class(tout) <- c("q2d_default", "data.frame")
  } else {
    if (justbest) tout <- tout$best
  }

  attr(tout, 'samp') <- samp
  attr(tout, 'criterion') <- criterion
  tout
}
