#' Normalize distribution name
#'
#' Maps various distribution name aliases to canonical short names
#' used by R's distribution functions.
#'
#' @param dist Character string naming a distribution.
#' @return Canonical distribution name (e.g. \code{"norm"}, \code{"lnorm"}).
#' @keywords internal
normalize_dist <- function(dist) {
  dist_map <- c(normal = 'norm', norm = 'norm', lnorm = 'lnorm', lnormal = 'lnorm',
                lognormal = 'lnorm', exp = 'exp', exponential = 'exp', beta = 'beta',
                gamma = 'gamma', weib = 'weibull', weibull = 'weibull')
  key <- tolower(dist)
  if(key %in% names(dist_map)) dist_map[[key]] else key
}

#' Get a distribution function by type
#'
#' Retrieves the d/p/q function for the named distribution, adding
#' a \code{...} formals entry so extra arguments are silently ignored.
#'
#' @param dist Character distribution name.
#' @param type One of \code{"d"}, \code{"p"}, \code{"q"}.
#' @return A function.
#' @keywords internal
get_dist_fn <- function(dist, type) {
  thisf <- get(paste0(type, normalize_dist(dist)))
  formals(thisf) <- c(formals(thisf), alist("..." = ))
  return(thisf)
}

#' @rdname get_dist_fn
#' @keywords internal
get_pfun <- function(dist) {
  get_dist_fn(dist, 'p')
}

#' @rdname get_dist_fn
#' @keywords internal
get_qfun <- function(dist) {
  get_dist_fn(dist, 'q')
}

#' @rdname get_dist_fn
#' @keywords internal
get_dfun <- function(dist) {
  get_dist_fn(dist, 'd')
}

# Ensure stats functions accessed via get() are declared as imports
#' @importFrom stats dnorm pnorm qnorm dlnorm plnorm qlnorm dexp pexp qexp
#'   dbeta pbeta qbeta dgamma pgamma qgamma dweibull pweibull qweibull
#'   integrate optim glm binomial pt residuals model.matrix
NULL
