#' Parse quantile input into internal data frame
#'
#' Converts a named numeric vector of quantiles (plus sample size \code{n})
#' into the internal data frame format used by the likelihood functions.
#' Handles quantile labels like \code{"25\%"}, \code{"50\%"}, \code{"min"},
#' \code{"max"}, as well as \code{"mean"}, \code{"sd"}, and \code{"se"}.
#'
#' @param samp Named numeric vector with quantile values and \code{n}.
#' @return A data frame with columns \code{n}, \code{lb}, \code{ub},
#'   \code{dtype}, \code{N}, \code{q}.
#'
#' @importFrom stats qnorm
#' @keywords internal
parse_quantiles <- function(samp) {
  thisn <- samp[['n']]
  n0 <- thisn-1

  thismeansd <- samp[(meansds <- which(tolower(names(samp)) %in% c("mean", "sd", "se")))]
  if(length(meansds)) samp <- samp[-meansds]


  xn <- names(samp)
  xn2 <- xn[-which(xn == "n")]
  thesequants <- samp[xn2]
  nquants <- length(thesequants)

  ordn0 <- suppressWarnings(as.numeric(gsub(pattern = "%", "", xn2))/100)*n0
  if("min" %in% xn2) ordn0[match(c("min", "max"), xn2)] <- c(0, n0)

  q <- ordn0/n0

  thisq <- thisn-nquants
  ordn <- ordn0+1

  tmp <- data.frame(
    n = c(rep(1, nquants), (c(ordn, thisn+1)-c(0, ordn))-1),
    lb =  c(thesequants, -Inf, thesequants),
    ub = c(thesequants, thesequants, Inf),
    dtype = rep(0:1, c(nquants, nquants+1)),
    N = rep(c(thisn, thisq), c(nquants, nquants+1)),
    q = c(rep(1/thisn, nquants), thisq * (c(q, 1) - c(0, q))))



  if(length(thismeansd)) {
    if("mean" %in% names(thismeansd)) {
      meansddf <- data.frame(

        n = thisn,
        lb = thismeansd[['mean']],
        ub = thismeansd[['mean']],
        dtype = 2,
        N = thisn,
        q = NA)
      if("sd" %in% names(thismeansd)) {
        meansddf$ub <- thismeansd[['sd']]
        meansddf$dtype <- 3
      }
      if("se" %in% names(thismeansd)) {
        meansddf$ub <- thismeansd[['se']]
        meansddf$dtype <- 4
      }

      tmp <- rbind(tmp, meansddf)
    }


  }
  tmp <- tmp[!(tmp$lb == -Inf & tmp$ub == Inf),]

  tmp

}
