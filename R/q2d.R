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
    weibull.scale.start = mean_hat / gamma((k_start + 1) / k_start)
  )
}

#' Estimate starting values for all distribution families
#'
#' Prepares starting parameter values for optimization by combining
#' \code{estmeansd::qe.fit} results with method-of-moments fallbacks.
#'
#' @param samp Named numeric vector with quantile values and \code{n}.
#' @param retorig Logical; if \code{TRUE}, return raw \code{qe.fit} output.
#' @param print_debug Logical; if \code{TRUE}, print debug information.
#' @return A named list of starting parameter vectors, one per distribution.
#' @keywords internal
estimate_start_values <- function(samp, retorig = FALSE, print_debug = FALSE) {
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
  list(normal = if(is.na(tmp$norm.par[1])) {
    c(mean = tmps$norm.mu.start[[1]], sd = tmps$norm.sigma.start[[1]])
  }
  else {
    c(mean = tmp$norm.par[['mu']], sd = tmp$norm.par[['sigma']])
  },
  lognormal = if(is.na(tmp$lnorm.par[1])) {
    c(meanlog = tmps$lnorm.mu.start[[1]], sdlog = tmps$lnorm.sigma.start[[1]])
  }
  else {
    c(meanlog = tmp$lnorm.par[['mu']], sdlog = tmp$lnorm.par[['sigma']])
  },
  exp = c(rate = 0.5),
  beta = if(is.na(tmp$beta.par[1])) {
    c(shape1 = tmps$beta.shape1.start[[1]], shape2 = tmps$beta.shape2.start[[1]])
  }
  else {
    c(tmp$beta.par)
  },
  scalebeta = if(is.na(tmp$beta.par[1])) {
    c(shape1 = tmps$beta.shape1.start[[1]], shape2 = tmps$beta.shape2.start[[1]], SCALE = 1)
  }
  else {
    c(tmp$beta.par, SCALE = 1)
  },
  gamma = if(is.na(tmp$gamma.par[1])) {
    c(rate = tmps$gamma.rate.start[[1]], shape = tmps$gamma.shape.start[[1]])
  }
  else {
    tmp$gamma.par
  },
  weibull = if(is.na(tmp$weibull.par[1])) {
    c(shape = tmps$weibull.shape.start[[1]], scale = tmps$weibull.scale.start[[1]])
  }
  else {
    tmp$weibull.par
  })

}


#' Fit parametric distributions to reported quantile data
#'
#' Fits six parametric distribution families (normal, lognormal, exponential,
#' beta, gamma, Weibull) to quantile-reported summary statistics using maximum
#' likelihood estimation, and selects the best-fitting distribution.
#'
#' @param samp A named numeric vector of quantile values. Names should be
#'   quantile labels such as \code{"25\%"}, \code{"50\%"}, \code{"min"},
#'   \code{"max"}, and must include \code{"n"} for sample size. May also
#'   include \code{"mean"}, \code{"sd"}, or \code{"se"}.
#' @param justbest Logical; if \code{TRUE} (default), return only the
#'   best-fitting distribution. If \code{FALSE}, return results for all six.
#' @param hessian Logical; if \code{TRUE}, compute the Hessian matrix at the
#'   optimum for each distribution.
#' @param tol Integer; number of decimal places for comparing log-likelihoods
#'   when selecting the best distribution (default 8).
#' @param return_df Logical; if \code{TRUE}, return results as a data frame.
#'
#' @return If \code{justbest = TRUE} and \code{return_df = FALSE} (default),
#'   a named numeric vector of moments, quantiles, and parameters for the
#'   best-fitting distribution. If \code{return_df = TRUE}, a data frame
#'   with one row per distribution (or just the best). The result has an
#'   attribute \code{"vector"} containing the original input.
#'
#' @examples
#' # Fit to IQR + median
#' q2d(c("25%" = 10, "50%" = 15, "75%" = 20, n = 100))
#'
#' # Fit to 5-number summary
#' q2d(c(min = 2, "25%" = 10, "50%" = 15, "75%" = 20, max = 35, n = 100))
#'
#' # Return all distributions as a data frame
#' q2d(c("25%" = 10, "50%" = 15, "75%" = 20, n = 100),
#'     justbest = FALSE, return_df = TRUE)
#'
#' @export
q2d <- function(samp, justbest = TRUE, hessian = FALSE, tol = 8, return_df = FALSE) {
  stvs <- estimate_start_values(samp)
  tmpout <- suppressWarnings(list(
    normal =        fit_dist(vals = samp, par = stvs$normal, dist = "normal", hessian = hessian),
    lognormal =     fit_dist(vals = samp, par = stvs$lognormal, dist = "lognormal", hessian = hessian),
    exponential =   fit_dist(vals = samp, par = stvs$exp, dist = "exp", hessian = hessian),
    beta =          fit_dist(vals = samp, par = stvs$beta, dist = "beta", hessian = hessian),
    gamma =         fit_dist(vals = samp, par = stvs$gamma, dist = "gamma", hessian = hessian),
    weibull =       fit_dist(vals = samp, par = stvs$weibull, dist = "weibull", hessian = hessian)
  ))

  tn <- names(tmpout)
  names(tn) <- tn

  tout <- lapply(tn, function(x) suppressWarnings(c(dist_moments(dist = x, pars = tmpout[[x]]$par, quantiles = c("0.1%" = 0.001, "2.5%" = 0.025, "5%" = 0.05, "25%" = 0.25, "50%" = 0.5, "75%" = 0.75, "95%" = 0.95, "97.5%" = 0.975, "99.9%" = 0.999)),
                                                    value = -tmpout[[x]]$value,
                                                    par1 = tmpout[[x]]$par[[1]],
                                                    par2 = unname(tmpout[[x]]$par[2]))))

  values <- sapply(tout, function(x) x[["value"]])
  minvs <- round(values-max(values, na.rm = TRUE), tol)
  minv <- match(0, minvs)
  tout$best <- c(tout[[minv]], "dist" = minv)
  tout$pars <- lapply(tn, function(x) c(tmpout[[x]]$par))
  tout$values <- values
  tout$n <- samp[['n']]
  if(return_df) {
    tmp <- suppressWarnings(as.data.frame(do.call(rbind, tout[1:6])))
    tmp$logLik <- tmp$value
    tmp$dist = rownames(tmp)
    tmp$best <- 0
    tmp$best[minv]<-1
    tout <- tmp[,c("dist", "par1", "par2", "logLik", "best",
                   "mean", "var", "sd", "median", "mode", "max",
                   "0.1%", "2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%", "99.9%")]
    if(justbest) tout <- tout[minv,]

  } else {
    if(justbest) tout <- tout$best
  }


  attr(tout, 'vector') <- samp
  tout

}
