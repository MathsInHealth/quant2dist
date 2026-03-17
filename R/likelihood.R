#' Joint probability for mean/SD observations
#'
#' Computes the joint log-probability of observing a mean (and optionally SD)
#' under the fitted distribution, using numerical integration.
#'
#' @param n Sample size.
#' @param dist1_pars Parameters for the first distribution.
#' @param dist2_pars Parameters for the second distribution.
#' @param dist1 Name of the first distribution.
#' @param dist2 Name of the second distribution.
#' @param log.p Logical; if \code{TRUE}, return log probability.
#' @return Numeric log-probability (or probability if \code{log.p = FALSE}).
#' @keywords internal
joint_prob <- function(n, dist1_pars, dist2_pars = c(mean = 0, sd = 1), dist1 = "normal", dist2 = "normal", log.p = TRUE) {

  if(inherits(tmp <- try(integrate(lower = 0, upper = 1, f =
                                function(x) do.call(get_dfun(dist1),
                                                    c(list(log = TRUE,
                                                           x = do.call(
                                                             get_qfun(dist2),
                                                             c(list(p = x),
                                                               as.list(dist2_pars)))),
                                                      as.list(dist1_pars)))), silent = TRUE), "try-error")) return(NA)

  tmp <- n*tmp[[1]]
  ifelse(log.p, tmp, exp(tmp))
}


#' Negative log-likelihood for quantile data
#'
#' Computes the negative log-likelihood of the data under a parametric
#' distribution, handling quantile intervals, discrete observations,
#' and reported mean/SD summaries.
#'
#' @param parnames Character vector of parameter names.
#' @param thesevals Data frame of parsed quantile data with parameter columns.
#' @param dist Character distribution name.
#' @return Scalar negative log-likelihood.
#'
#' @importFrom stats dnorm
#' @keywords internal
neg_log_lik <- function(parnames, thesevals, dist = "normal") {
  pfun <- get_pfun(dist)
  dfun <- get_dfun(dist)

  -(sum(    do.call(dfun, c(list(thesevals[thesevals$dtype==0,"lb"], log   = TRUE),  as.list(unname(thesevals[thesevals$dtype==0,parnames, drop = FALSE])))) *   # Discrete observations
              thesevals[thesevals$dtype==0,"n"])
    +                                                               # multiplied by n
      sum(log(pmax(do.call(pfun, c(list(thesevals[thesevals$dtype==1,"ub"], log.p = FALSE),  as.list(unname(thesevals[thesevals$dtype==1,parnames, drop = FALSE]))))    - # CDF up to upper bound
                     do.call(pfun, c(list(thesevals[thesevals$dtype==1,"lb"], log.p = FALSE),  as.list(unname(thesevals[thesevals$dtype==1,parnames, drop = FALSE])))),
                   rep(.Machine$double.xmin, sum(thesevals$dtype == 1)))) * # minus CDF up to lower bound
            thesevals[thesevals$dtype == 1,"n"]) +                        # multiplied by n
      sum(unlist(lapply(which(thesevals$dtype == 2), function(thisv) {    # Reported mean, just one case
        do.call(dnorm, c(list(thesevals$lb[thisv], log = TRUE), mean_se(dist = dist, par = as.list(unname(thesevals[thisv, parnames, drop = FALSE])), n = thesevals[thisv, "n"])))
      }))) +
      sum(unlist(lapply(which(thesevals$dtype == 3), function(thisv) {    # Reported mean + SD, n cases
        joint_prob(n = thesevals$n[thisv], dist1_pars = mean_sd(dist = dist, par = unname(thesevals[thisv, parnames, drop = FALSE])), dist2_pars = c(mean = thesevals[thisv, "lb"], sd = thesevals[thisv, "ub"]))
      }))) +
      sum(unlist(lapply(which(thesevals$dtype == 4), function(thisv) {    # Reported mean + SE, just one case
        joint_prob(n = 1, dist1_pars = mean_se(dist = dist, par = unname(thesevals[thisv, parnames, drop = FALSE]), n = thesevals[thisv, "n"]), dist2_pars = c(mean = thesevals[thisv, "lb"], sd = thesevals[thisv, "ub"]))
      }))) +
      sum(unlist(lapply(which(thesevals$dtype == 5), function(thisv) {    # Reported arguments for dist
        joint_prob(n = thesevals$n[thisv], dist1 = dist, dist2 = dist, dist1_pars = unname(thesevals[thisv, parnames, drop = FALSE]), dist2_pars = c(thesevals[thisv, "lb"], thesevals[thisv, "ub"]))
      })))
  )
}

#' Wrapper for neg_log_lik used by optimizer
#'
#' Binds parameter values into the data frame and calls \code{neg_log_lik}.
#'
#' @param par Named numeric vector of distribution parameters.
#' @param thesevals Data frame of parsed quantile data.
#' @param dist Character distribution name.
#' @return Scalar negative log-likelihood.
#' @keywords internal
nll_wrap <- function(par, thesevals, dist) {
  neg_log_lik(names(par), cbind(thesevals, as.list(par)), dist)
}


#' Fit a single distribution to quantile data
#'
#' Parses the input, then minimizes the negative log-likelihood using
#' \code{\link[ucminf]{ucminf}}.
#'
#' @param vals Named numeric vector of quantiles (see \code{\link{q2d}}).
#' @param par Named numeric starting values for the distribution parameters.
#' @param dist Character distribution name.
#' @param hessian Logical; compute Hessian at the optimum?
#' @param ret_data Logical; if \code{TRUE}, return parsed data instead of fitting.
#' @return An object returned by \code{\link[ucminf]{ucminf}}.
#'
#' @importFrom ucminf ucminf
#' @keywords internal
fit_dist <- function(vals, par = c(mean = 0, sd = 1), dist = "normal", hessian = TRUE, ret_data = FALSE) {

  thesevals <- parse_quantiles(vals)

  if(ret_data) return(thesevals)

  ucminf(par = par, fn = nll_wrap, thesevals = thesevals, dist = dist, hessian = hessian)

}
