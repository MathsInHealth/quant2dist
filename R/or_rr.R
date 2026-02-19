#' Compute negative Poisson log-likelihood via integration
#'
#' Helper function for \code{\link{estimate_rr}} that computes the negative
#' log-likelihood for the robust Poisson (log-link) regression model by
#' integrating over the covariate distribution.
#'
#' @param par Named numeric vector with elements \code{INTERCEPT} and \code{QOI}
#'   (the log relative risk per unit increase in X).
#' @param distdf Data frame containing distribution parameters (see \code{\link{estimate_rr}}).
#'
#' @return Numeric scalar: the negative log-likelihood.
#'
#' @details
#' For each group in \code{distdf}, the log-likelihood contribution is computed by
#' integrating \code{y*log(mu) - mu} (where \code{mu = exp(intercept + beta*X)})
#' over the quantile function of X, then scaling by sample size n.
#'
#' @keywords internal
#' @importFrom stats integrate
rr_log_likelihood <- function(par, distdf) {
  outv <- -do.call(sum, (lapply(1:NROW(distdf), function(x) {
    qfun <- get_qfun(tolower(distdf$dist[x]))
    tx <- try(integrate(lower = 0, upper = 1, f = function(y) {
      
      Xval <- do.call(qfun, c(list(p = y), unname(as.list(distdf[x, c("par1", "par2")]))))
      Xb <- par['INTERCEPT'] + par['QOI'] * Xval
      
      mu <- exp(Xb)
      
      # Poisson log-likelihood (working): y*log(mu) - mu
      # Clamp mu to avoid log(0) or overflow
      mu <- pmax(mu, .Machine$double.xmin)
      mu <- pmin(mu, .Machine$double.xmax / 2)
      
      ll_val <- distdf$outcome[x] * log(mu) - mu
      
      ll_val[!is.finite(ll_val)] <- log(.Machine$double.xmin)
      ll_val
    }, subdivisions = 200L), silent = TRUE)
    
    if (inherits(tx, 'try-error')) {
      log(.Machine$double.xmin)
    } else {
      distdf$n[x] * tx[[1]]
    }
  })))
  outv
}

#' Compute sandwich variance meat matrix via integration
#'
#' Helper function for \code{\link{estimate_rr}} that computes the "meat"
#' component of the sandwich variance estimator by integrating score outer
#' products over covariate distributions.
#'
#' @param par Named numeric vector with elements \code{INTERCEPT} and \code{QOI}.
#' @param distdf Data frame containing distribution parameters (see \code{\link{estimate_rr}}).
#'
#' @return Matrix (2x2): the meat component \code{B} in the sandwich formula
#'   \code{V = A^{-1} B A^{-1}}.
#'
#' @details
#' For the robust Poisson model, the meat matrix is:
#' \code{B = sum_i n_i * E[ (y_i - mu(X))^2 * z z' ]}
#' where \code{z = (1, X)'} and expectations are computed via integration
#' over the quantile function of X.
#'
#' @keywords internal
#' @importFrom stats integrate
rr_sandwich_meat <- function(par, distdf) {
  p <- length(par)
  meat <- matrix(0, p, p)
  
  for (x in 1:NROW(distdf)) {
    qfun <- get_qfun(tolower(distdf$dist[x]))
    
    # We need to integrate each element of the 2x2 outer product matrix.
    # Since p=2, we integrate 3 unique elements: (1,1), (1,2), (2,2).
    
    meat_elements <- sapply(1:p, function(j) {
      sapply(1:p, function(k) {
        tx <- try(integrate(lower = 0, upper = 1, f = function(y) {
          Xval <- do.call(qfun, c(list(p = y), unname(as.list(distdf[x, c("par1", "par2")]))))
          Xb <- par['INTERCEPT'] + par['QOI'] * Xval
          
          mu <- exp(Xb)
          mu <- pmax(mu, .Machine$double.xmin)
          mu <- pmin(mu, .Machine$double.xmax / 2)
          
          resid <- distdf$outcome[x] - mu
          
          # z = (1, Xval)
          z_j <- if (j == 1) 1 else Xval
          z_k <- if (k == 1) 1 else Xval
          
          resid^2 * z_j * z_k
        }, subdivisions = 200L), silent = TRUE)
        
        if (inherits(tx, 'try-error')) 0 else tx[[1]]
      })
    })
    
    # meat_elements is p x p
    # Scale by n_i^2 because within a group of n_i identical-outcome obs,
    # the summed score is n_i times the per-obs score, and the outer 
    # product of the summed score is n_i^2 times the per-obs outer product.
    # But we then have n_i such "clusters" ... Actually, in this grouped
    # setup each row IS the cluster. The total score contribution from
    # group i is n_i * E[score_per_obs], and its outer product contributes
    # n_i^2 * E[score_per_obs %o% score_per_obs] to the meat.
    #
    # However, the standard robust Poisson (as in individual-level data)
    # sums over *individual* score outer products. With n_i individuals
    # in group i, each having the same score integrand, this gives:
    # n_i * E[(y - mu)^2 * zz']
    #
    # The distinction matters: are we treating this as n_i independent 
    # observations (each with their own X drawn from the distribution),
    # or as one "cluster" of n_i identical observations?
    #
    # In the simulation approach (exp_OR_RR), each simulated observation
    # is treated independently. The robust Poisson GLM on individual data
    # uses HC (heteroskedasticity-consistent) SEs, NOT clustered SEs.
    # So each obs contributes its own score outer product.
    #
    # Therefore the correct scaling is n_i * integral (not n_i^2):
    
    meat <- meat + distdf$n[x] * meat_elements
  }
  
  meat
}

#' Estimate relative risks from distribution parameters via integration
#'
#' Computes relative risks (RR) relating a continuous predictor to a binary
#' outcome using the robust Poisson method (log-link GLM with sandwich-corrected
#' standard errors). Unlike simulation-based approaches, this function integrates
#' the log-likelihood directly over the predictor distributions, providing
#' more precise and faster estimates.
#'
#' @param distdf A data frame with one row per outcome group, containing:
#'   \itemize{
#'     \item \code{outcome}: Binary outcome indicator (0 or 1).
#'     \item \code{dist}: Distribution name (e.g., \code{"normal"}, \code{"lognormal"}).
#'     \item \code{par1}: First distribution parameter (e.g., mean or shape).
#'     \item \code{par2}: Second distribution parameter (e.g., SD or scale).
#'     \item \code{n}: Sample size for each group.
#'     \item \code{mean}, \code{sd}: Optional; if not provided, computed from \code{par1}/\code{par2}.
#'   }
#' @param colnms Named character vector mapping standard column names to
#'   actual column names in \code{distdf}. Default assumes standard names.
#' @param stvs Optional starting values for optimization: \code{c(intercept, slope)}.
#' @param ll_fun Log-likelihood function (default: \code{\link{rr_log_likelihood}}).
#'
#' @return A named numeric vector with components:
#'   \item{RR}{Relative risk (exponentiated slope coefficient).}
#'   \item{log(RR)}{Log relative risk (slope coefficient).}
#'   \item{SE(log(RR))}{Sandwich-adjusted standard error of log(RR).}
#'   \item{p}{Two-sided p-value testing log(RR) = 0.}
#'   \item{LB2.5, UB97.5}{95\% confidence interval for RR.}
#'
#'   The vector has an attribute \code{"complete"} containing the full
#'   \code{ucminf} output, including \code{vcov} (robust sandwich variance),
#'   \code{vcov_naive} (model-based variance), and \code{SE_naive}.
#'
#' @details
#' This function uses the robust Poisson regression approach, which models
#' \code{P(Y=1|X) = exp(alpha + beta*X)} using a Poisson quasi-likelihood.
#' The sandwich variance estimator corrects for variance misspecification.
#'
#' The method integrates the log-likelihood over the quantile functions of X
#' in each outcome group, avoiding Monte Carlo error from simulation-based methods.
#'
#' Supported distributions: normal, lognormal, exponential, beta, gamma, weibull.
#'
#' @examples
#' \dontrun{
#' # Example: RR for normal distributions
#' df <- data.frame(
#'   outcome = c(0, 1, 0, 1),
#'   dist = "normal",
#'   par1 = c(9.3, 8.7, 8.4, 7.7),  # means
#'   par2 = c(1.6, 1.9, 1.3, 1.3),  # SDs
#'   n = c(188, 38, 188, 38)
#' )
#' result <- estimate_rr(df)
#' print(result)
#'
#' # Access robust vs naive SEs
#' attr(result, "complete")$SE         # sandwich-adjusted
#' attr(result, "complete")$SE_naive   # model-based
#' }
#'
#' @seealso \code{\link{estimate_or}} for odds ratios
#' @importFrom stats qnorm pt
#' @importFrom ucminf ucminf
#' @export
estimate_rr <- function(distdf,
                   colnms = c(outcome = 'outcome', dist = 'dist', par1 = 'par1', 
                              par2 = 'par2', n = 'n', mean = 'mean', sd = 'sd'),
                   stvs = NULL,
                   ll_fun = rr_log_likelihood) {
  
  thesecolnms <- c(outcome = 'outcome', dist = 'dist', par1 = 'par1', 
                   par2 = 'par2', n = 'n', mean = 'mean', sd = 'sd')
  colnms <- colnms[names(colnms) %in% names(thesecolnms)]
  if (length(colnms)) {
    thesecolnms[names(colnms)] <- colnms
  }
  
  if (!all(thesecolnms[1:5] %in% colnames(distdf))) {
    stop('distdf does not contain necessary information.')
  }
  
  if (!all(thesecolnms[6:7] %in% colnames(distdf))) {
    distdf[, thesecolnms[6:7]] <- t(sapply(1:NROW(distdf), function(x) {
      moments <- dist_moments(dist = distdf[x, thesecolnms[['dist']]],
                              pars = unname(as.vector(as.matrix(distdf[x, thesecolnms[c('par1', "par2")]]))),
                              included = c('mean', 'sd'),
                              quantiles = NULL)
      moments[c('mean', 'sd')]
    }))
  }
  distdf <- distdf[, thesecolnms]
  colnames(distdf) <- names(thesecolnms)
  
  nsum <- sum(distdf$n)
  
  # Starting values for log link:
  # Intercept: log of overall event rate
  event_rate <- sum(distdf$n[distdf$outcome == 1]) / nsum
  icpt <- log(max(event_rate, .Machine$double.xmin))
  QOI <- 0
  
  if (!is.null(stvs)) {
    icpt <- stvs[[1]]
    QOI <- stvs[[2]]
  }
  
  # Optimize Poisson (working) log-likelihood
  oout <- ucminf(par = c(INTERCEPT = icpt, QOI = QOI),
                 fn = ll_fun,
                 distdf = distdf,
                 hessian = 1)
  
  oout$SE <- rep(NA, length(oout$par))
  names(oout$SE) <- names(oout$par)
  
  # --- Model-based (naive) variance from Hessian ---
  bread <- solve(oout$hessian)  # A^{-1}
  
  # --- Robust sandwich variance ---
  meat <- rr_sandwich_meat(oout$par, distdf)
  
  oout$vcov_naive <- bread
  oout$vcov <- bread %*% meat %*% bread   # sandwich
  oout$SE[] <- sqrt(diag(oout$vcov))
  
  oout$SE_naive <- sqrt(diag(oout$vcov_naive))
  
  oout$t <- abs(oout$par) / oout$SE
  oout$p <- 2 * pt(oout$t, df = nsum - 1, lower.tail = FALSE)
  
  # Output
  tmp <- c(RR             = exp(oout$par[[2]]),
           "log(RR)"      = oout$par[[2]],
           "SE(log(RR))"  = oout$SE[[2]],
           p              = oout$p[[2]],
           'LB2.5'        = exp(oout$par[[2]] + qnorm(0.025) * oout$SE[[2]]),
           'UB97.5'       = exp(oout$par[[2]] + qnorm(0.975) * oout$SE[[2]]))
  attr(tmp, 'complete') <- oout
  tmp
}

#' Compute negative logistic log-likelihood via integration
#'
#' Helper function for \code{\link{estimate_or}} that computes the negative
#' log-likelihood for logistic regression by integrating over the covariate
#' distribution.
#'
#' @param par Named numeric vector with elements \code{INTERCEPT} and \code{QOI}
#'   (the log odds ratio per unit increase in X).
#' @param distdf Data frame containing distribution parameters (see \code{\link{estimate_or}}).
#'
#' @return Numeric scalar: the negative log-likelihood.
#'
#' @details
#' For each group in \code{distdf}, the log-likelihood contribution is computed by
#' integrating the logistic log-likelihood over the quantile function of X,
#' then scaling by sample size n. Uses \code{tanh} for numerical stability.
#'
#' @keywords internal
#' @importFrom stats integrate
or_log_likelihood <- function(par, distdf) {
  outv <- -do.call(sum, (lapply(1:NROW(distdf), function(x) {
    qfun <- get_qfun(tolower(distdf$dist[x]))
    tx <- try(integrate(lower = 0, upper = 1, f = function(y) {
      
      Xb <- (par['INTERCEPT']+par['QOI']*do.call(qfun, c(list(p = y), unname(as.list(distdf[x, c("par1", "par2")])))))
      
      logistic_tmp <- .5+.5*tanh(Xb/2)
      pvals <- log(distdf$outcome[x] *    logistic_tmp +
                     (1-distdf$outcome[x])* (1-logistic_tmp))
      
      pvals[pvals == -Inf] <- log(.Machine$double.xmin)
      pvals
    }, subdivisions = 200L), silent = TRUE)
    if(inherits(tx, 'try-error')) {
      log(.Machine$double.xmin)
    } else {
      distdf$n[x]*tx[[1]]
    }
  })))
  outv
}

#' Estimate odds ratios from distribution parameters via integration
#'
#' Computes odds ratios (OR) relating a continuous predictor to a binary
#' outcome using logistic regression with log-likelihood integrated directly
#' over the predictor distributions. This integration-based approach is more
#' precise than simulation methods and avoids Monte Carlo error.
#'
#' @param distdf A data frame with one row per outcome group, containing:
#'   \itemize{
#'     \item \code{outcome}: Binary outcome indicator (0 or 1).
#'     \item \code{dist}: Distribution name (e.g., \code{"normal"}, \code{"lognormal"}).
#'     \item \code{par1}: First distribution parameter (e.g., mean or shape).
#'     \item \code{par2}: Second distribution parameter (e.g., SD or scale).
#'     \item \code{n}: Sample size for each group.
#'     \item \code{mean}, \code{sd}: Optional; if not provided, computed from \code{par1}/\code{par2}.
#'   }
#' @param colnms Named character vector mapping standard column names to
#'   actual column names in \code{distdf}. Default assumes standard names.
#' @param stvs Optional starting values for optimization: \code{c(intercept, slope)}.
#' @param ll_fun Log-likelihood function (default: \code{\link{or_log_likelihood}}).
#'
#' @return A named numeric vector with components:
#'   \item{OR}{Odds ratio (exponentiated slope coefficient).}
#'   \item{log(OR)}{Log odds ratio (slope coefficient).}
#'   \item{SE(log(OR))}{Standard error of log(OR) from model-based variance.}
#'   \item{p}{Two-sided p-value testing log(OR) = 0.}
#'   \item{LB2.5, UB97.5}{95\% confidence interval for OR.}
#'
#'   The vector has an attribute \code{"complete"} containing the full
#'   \code{ucminf} output, including \code{vcov} (variance-covariance matrix)
#'   and model parameters.
#'
#' @details
#' This function fits a logistic regression model by integrating the log-likelihood
#' over the quantile functions of the predictor X in each outcome group. The
#' integration is performed using adaptive quadrature (\code{stats::integrate}).
#'
#' Starting values are computed from the log odds of outcome prevalence and
#' a standardized mean difference-based approximation.
#'
#' Supported distributions: normal, lognormal, exponential, beta, gamma, weibull.
#'
#' @examples
#' \dontrun{
#' # Example: OR for normal distributions
#' df <- data.frame(
#'   outcome = c(0, 1, 0, 1),
#'   dist = "normal",
#'   par1 = c(9.3, 8.7, 8.4, 7.7),  # means
#'   par2 = c(1.6, 1.9, 1.3, 1.3),  # SDs
#'   n = c(188, 38, 188, 38)
#' )
#' result <- estimate_or(df)
#' print(result)
#' }
#'
#' @seealso \code{\link{estimate_rr}} for relative risks
#' @importFrom stats qnorm pt
#' @importFrom ucminf ucminf
#' @export
estimate_or <- function(distdf,
                   colnms = c(outcome = 'outcome', dist = 'dist', par1 = 'par1', par2 = 'par2', n = 'n', mean = 'mean', sd = 'sd'),
                   stvs = NULL,
                   ll_fun = or_log_likelihood) {
  
  thesecolnms <- c(outcome = 'outcome', dist = 'dist', par1 = 'par1', par2 = 'par2', n = 'n', mean = 'mean', sd = 'sd')
  colnms <- colnms[names(colnms) %in% names(thesecolnms)]
  if(length(colnms)) {
    thesecolnms[names(colnms)] <- colnms
  }
  if(!all(thesecolnms[1:5] %in% colnames(distdf))) stop('distdf does not contain necessary information.')
  
  if(!all(thesecolnms[6:7] %in% colnames(distdf))) { # Missing mean/sd
    distdf[, thesecolnms[6:7]] <- t(sapply(1:NROW(distdf), function(x) {
          moments <- dist_moments(dist = distdf[x, thesecolnms[['dist']]],
                                  pars = unname(as.vector(as.matrix(distdf[x, thesecolnms[c('par1', "par2")]]))),
                                  included = c('mean', 'sd'),
                                  quantiles = NULL)
          moments[c('mean', 'sd')]
    }))
  }
  distdf <- distdf[, thesecolnms]
  colnames(distdf) <- names(thesecolnms)
  
  nsum <- sum(distdf$n)
  
  smd <-((mean(distdf$mean[distdf$outcome == 1])-mean(distdf$mean[distdf$outcome == 0]))/
           sqrt(sum(distdf$n * distdf$sd^2)/
                  (sum(distdf$n))))
  
  QOI <- (max(min(smd*pi/sqrt(3), 0.5), -0.5))
  QOI <- 0
  
  icpt <- log(sum(distdf$n[distdf$outcome == 1])/
                sum(distdf$n[distdf$outcome == 0]))
  
  if(!is.null(stvs)) {
    icpt <- stvs[[1]]
    QOI <- stvs[[2]]
  }
  
  oout <- ucminf(par = c(INTERCEPT = icpt,
                         QOI = QOI),
                 fn = ll_fun,
                 distdf = distdf,
                 hessian = 1)

  oout$SE <- rep(NA, length(oout$par))
  names(oout$SE) <- names(oout$par)
  
  
  
  
  oout$vcov <-     solve(oout$hessian)
  oout$SE[] <- sqrt(diag(oout$vcov))
  
  oout$t <- abs(oout$par)/oout$SE
  oout$p <- 2 * pt(oout$t, df = nsum - 1, lower.tail = FALSE)
  
  
  
  tmp <-c(OR = exp(oout$par[[2]]),
          "log(OR)" = oout$par[[2]],
          "SE(log(OR))" = oout$SE[[2]],
          p = oout$p[[2]],
          'LB2.5' = exp(oout$par[[2]]+qnorm(0.025)*oout$SE[[2]]),
          'UB97.5' = exp(oout$par[[2]]+qnorm(0.975)*oout$SE[[2]]))
  attr(tmp, 'complete') <- oout
  tmp

}

#' Estimate OR and RR from reported quantiles
#'
#' Convenience wrapper that fits distributions to reported quantiles using
#' \code{\link{q2d}}, then estimates odds ratios and relative risks using
#' \code{\link{estimate_or}} and \code{\link{estimate_rr}}. This combines
#' the two-step workflow into a single function call.
#'
#' @param case_quants Named numeric vector of quantiles for cases (outcome=1),
#'   or a list of such vectors for multiple case groups/studies.
#'   Names should be quantile labels (e.g., \code{"min"}, \code{"25\%"},
#'   \code{"50\%"}, \code{"75\%"}, \code{"max"}) plus \code{"n"} for sample size.
#'   See \code{\link{q2d}} for supported formats.
#' @param control_quants Named numeric vector of quantiles for controls (outcome=0),
#'   or a list of such vectors for multiple control groups/studies.
#'   Same format as \code{case_quants}.
#' @param measure Character string: \code{"OR"} for odds ratio only, \code{"RR"}
#'   for relative risk only, or \code{"both"} (default) for both measures.
#' @param dist_family Optional: force a specific distribution family instead of
#'   automatic selection. Must be one of: \code{"normal"}, \code{"lognormal"},
#'   \code{"exponential"}, \code{"beta"}, \code{"gamma"}, \code{"weibull"}.
#'   If \code{NULL} (default), best-fitting distribution is selected for each group.
#'
#' @return A list with components:
#'   \item{case_fit}{Distribution fit(s) for cases from \code{\link{q2d}}.
#'     A single fit object if input is a vector, or a list of fit objects if input is a list.}
#'   \item{control_fit}{Distribution fit(s) for controls from \code{\link{q2d}}.
#'     A single fit object if input is a vector, or a list of fit objects if input is a list.}
#'   \item{distdf}{Data frame passed to OR/RR estimation functions.}
#'   \item{OR}{Odds ratio results (if computed), from \code{\link{estimate_or}}.}
#'   \item{RR}{Relative risk results (if computed), from \code{\link{estimate_rr}}.}
#'
#' @details
#' This function streamlines the workflow of:
#' \enumerate{
#'   \item Fitting parametric distributions to reported quantiles for cases and controls
#'   \item Extracting fitted parameters
#'   \item Computing OR and/or RR via integration-based methods
#' }
#'
#' If \code{dist_family} is not specified, the function selects the best-fitting
#' distribution separately for cases and controls. This allows for different
#' distributions in each group (e.g., normal for controls, lognormal for cases).
#'
#' When inputs are lists (multiple studies/groups), each element is fitted separately
#' and all groups are pooled in the final OR/RR estimation. This allows combining
#' data from multiple sources with potentially different distributions.
#'
#' @examples
#' \dontrun{
#' # Example 1: Single case and control group
#' cases <- c("min" = 2.1, "25%" = 4.5, "50%" = 6.2,
#'            "75%" = 8.3, "max" = 15.2, "n" = 38)
#' controls <- c("min" = 3.2, "25%" = 5.1, "50%" = 7.1,
#'               "75%" = 9.1, "max" = 14.5, "n" = 188)
#'
#' result <- q2d_or_rr(cases, controls)
#' print(result$OR)
#' print(result$RR)
#'
#' # Example 2: Multiple case groups (pooled analysis)
#' cases_list <- list(
#'   c("min" = 2.1, "25%" = 4.5, "50%" = 6.2, "75%" = 8.3, "max" = 15.2, "n" = 38),
#'   c("25%" = 3.7, "50%" = 6.0, "75%" = 7.5, "n" = 150)
#' )
#' controls <- c("min" = 3.2, "25%" = 5.1, "50%" = 7.1,
#'               "75%" = 9.1, "max" = 14.5, "n" = 188)
#'
#' result_pooled <- q2d_or_rr(cases_list, controls)
#' # distdf will have 3 rows: 1 control group + 2 case groups
#' print(result_pooled$distdf)
#' }
#'
#' @seealso \code{\link{q2d}}, \code{\link{estimate_or}}, \code{\link{estimate_rr}}
#' @export
q2d_or_rr <- function(case_quants, control_quants, measure = "both",
                      dist_family = NULL) {

  # Helper function to validate a single quantile vector
  validate_quants <- function(quants, name) {
    if (!all(c("n") %in% names(quants))) {
      stop(name, " must include 'n' (sample size)")
    }
  }

  # Check if inputs are lists or single vectors
  case_is_list <- is.list(case_quants) && !is.null(names(case_quants[[1]]))
  control_is_list <- is.list(control_quants) && !is.null(names(control_quants[[1]]))

  # Validate inputs
  if (case_is_list) {
    lapply(seq_along(case_quants), function(i) {
      validate_quants(case_quants[[i]], paste0("case_quants[[", i, "]]"))
    })
  } else {
    validate_quants(case_quants, "case_quants")
  }

  if (control_is_list) {
    lapply(seq_along(control_quants), function(i) {
      validate_quants(control_quants[[i]], paste0("control_quants[[", i, "]]"))
    })
  } else {
    validate_quants(control_quants, "control_quants")
  }

  measure <- tolower(measure)
  if (!measure %in% c("or", "rr", "both")) {
    stop("measure must be 'OR', 'RR', or 'both'")
  }

  # Fit distributions to quantiles
  if (case_is_list) {
    case_fit <- lapply(case_quants, function(q) q2d(q, justbest = FALSE))
  } else {
    case_fit <- q2d(case_quants, justbest = FALSE)
  }

  if (control_is_list) {
    control_fit <- lapply(control_quants, function(q) q2d(q, justbest = FALSE))
  } else {
    control_fit <- q2d(control_quants, justbest = FALSE)
  }

  # Helper function to extract parameters from a single fit
  extract_params <- function(fit, quants, dist_family = NULL) {
    if (is.null(dist_family)) {
      # Find best fit distribution (highest log-likelihood)
      dist <- names(which.max(fit$values))
      params <- fit[[dist]]
    } else {
      # Use specified distribution
      dist_key <- tolower(dist_family)
      if (!dist_key %in% names(fit)) {
        stop("Distribution '", dist_family, "' not available in fits. Available: ",
             paste(setdiff(names(fit), c("best", "pars", "values", "n")), collapse=", "))
      }
      dist <- dist_key
      params <- fit[[dist_key]]
    }

    list(
      dist = normalize_dist(dist),
      par1 = params["par1"],
      par2 = params["par2"],
      n = quants[["n"]]
    )
  }

  # Build distdf for OR/RR estimation
  # Control groups (outcome = 0)
  if (control_is_list) {
    control_rows <- lapply(seq_along(control_fit), function(i) {
      p <- extract_params(control_fit[[i]], control_quants[[i]], dist_family)
      data.frame(
        outcome = 0,
        dist = p$dist,
        par1 = unname(p$par1),
        par2 = unname(p$par2),
        n = unname(p$n),
        stringsAsFactors = FALSE
      )
    })
    control_df <- do.call(rbind, control_rows)
  } else {
    control_params <- extract_params(control_fit, control_quants, dist_family)
    control_df <- data.frame(
      outcome = 0,
      dist = control_params$dist,
      par1 = unname(control_params$par1),
      par2 = unname(control_params$par2),
      n = unname(control_params$n),
      stringsAsFactors = FALSE
    )
  }

  # Case groups (outcome = 1)
  if (case_is_list) {
    case_rows <- lapply(seq_along(case_fit), function(i) {
      p <- extract_params(case_fit[[i]], case_quants[[i]], dist_family)
      data.frame(
        outcome = 1,
        dist = p$dist,
        par1 = unname(p$par1),
        par2 = unname(p$par2),
        n = unname(p$n),
        stringsAsFactors = FALSE
      )
    })
    case_df <- do.call(rbind, case_rows)
  } else {
    case_params <- extract_params(case_fit, case_quants, dist_family)
    case_df <- data.frame(
      outcome = 1,
      dist = case_params$dist,
      par1 = unname(case_params$par1),
      par2 = unname(case_params$par2),
      n = unname(case_params$n),
      stringsAsFactors = FALSE
    )
  }

  # Combine control and case rows
  distdf <- rbind(control_df, case_df)
  rownames(distdf) <- NULL

  # Estimate OR and/or RR
  result <- list(
    case_fit = case_fit,
    control_fit = control_fit,
    distdf = distdf
  )

  if (measure %in% c("or", "both")) {
    result$OR <- estimate_or(distdf)
  }

  if (measure %in% c("rr", "both")) {
    result$RR <- estimate_rr(distdf)
  }

  result
}

# Test code removed - moved to tests/testthat/test-or_rr.R
