#' Compute mean and SD for a distribution
#'
#' Analytically computes the mean and standard deviation from distribution
#' parameters using closed-form expressions for each supported family.
#'
#' @param dist Character distribution name.
#' @param par Named numeric vector of distribution parameters.
#' @return A list with elements \code{mean} and \code{sd}.
#' @keywords internal
mean_sd <- function(dist = "normal", par = NULL) {
  dist <- normalize_dist(dist)
  if(is.null(par)) return(c(mean = NA, SD = NA))
  if(is.null(names(par))) {
    names(par) <- switch(dist,
                         normal =      c("mean", "sd"),
                         norm =        c("mean", "sd"),
                         lognormal =   c("meanlog", "sdlog"),
                         lnorm =       c("meanlog", "sdlog"),
                         exp =         "rate",
                         exponential = "rate",
                         beta =        c("shape1", "shape2"),
                         gamma =       c("shape", "rate"),
                         weibull =     c("shape", "scale"))

  }
  if(dist == "scalebeta" && !"SCALE" %in% names(par)) par$SCALE <- 1
  par <- as.list(par)
  if(any(is.na(par))) dist <- ""
  if(dist %in% c("exp", "exponential")) if(par$rate <=0) dist <- ""
  switch(dist,
         normal = list(mean = par$mean,
                       sd = par$sd),
         norm =   list(mean = par$mean,
                       sd = par$sd),
         lognormal = list(mean = exp(par$mean + 0.5* par$sd^2),
                          sd = sqrt(exp(2*(par$mean + par$sd^2))-exp(2*par$mean+par$sd^2))),
         lnormal = list(mean = exp(par$mean + 0.5* par$sd^2),
                        sd = sqrt(exp(2*(par$mean + par$sd^2))-exp(2*par$mean+par$sd^2))),
         lnorm = list(mean = exp(par$mean + 0.5* par$sd^2),
                      sd = sqrt(exp(2*(par$mean + par$sd^2))-exp(2*par$mean+par$sd^2))),
         exponential = list(mean = 1/par$rate,
                            sd = sqrt(1/(par$rate)^2)),
         exp = list(mean = 1/par$rate,
                    sd = sqrt(1/(par$rate)^2)),
         beta = list(mean = (par$shape1/(par$shape1+par$shape2)),
                     sd = sqrt((par$shape1*par$shape2)/(((par$shape1 + par$shape2)^2)*(par$shape1 + par$shape2 + 1)))),
         scalebeta = list(mean = (par$shape1/(par$shape1+par$shape2))/par$SCALE,
                          sd = sqrt((par$shape1*par$shape2)/(((par$shape1 + par$shape2)^2)*(par$shape1 + par$shape2 + 1))/par$SCALE)),
         gamma = list(mean = par$shape/par$rate,
                      sd = sqrt(par$shape/par$rate^2)),
         weibull = list(mean = par$scale*gamma(1+1/par$shape),
                        sd = sqrt((par$scale^2)*(gamma(1+2/par$shape)-(gamma(1+1/par$shape))^2))),
         list(mean = NA, sd = NA)
  )
}

#' Compute mean and SE for a distribution
#'
#' Computes the mean and standard error (SD / sqrt(n)) from distribution
#' parameters.
#'
#' @param dist Character distribution name.
#' @param par Named numeric vector of distribution parameters.
#' @param n Sample size.
#' @return A list with elements \code{mean} and \code{sd} (the SE).
#' @keywords internal
mean_se <- function(dist = "normal", par = NULL, n) {
  tmp <- mean_sd(dist, par)
  list(mean = tmp$mean, sd = tmp$sd/sqrt(n))
}


#' Compute distribution moments and quantiles
#'
#' Calculates mean, variance, SD, median, mode, and selected quantiles
#' from distribution parameters.
#'
#' @param dist Character distribution name.
#' @param pars Numeric vector of distribution parameters.
#' @param included Character vector of moment names to include.
#' @param quantiles Named numeric vector of quantile probabilities.
#' @return Named numeric vector of moments and quantiles.
#' @keywords internal
dist_moments <- function(dist = "normal", pars = NULL, included = c("mean", "var", "sd", "median", "mode"), quantiles = c("2.5%" = 0.025, "5%" = 0.05, "25%" = 0.25, "50%" = 0.5, "75%" = 0.75, "95%" = 0.95, "97.5%" = 0.975)) {
  dist <- tolower(dist)
  if(is.null(pars)) return(c(mean = NA, SD = NA))

  dist <- normalize_dist(dist)

  pars <- as.list(pars)
  names(pars) <- paste0('par', seq_along(pars))
  if(!is.null(quantiles)) {
    names(quantiles) <- paste0(quantiles*100, "%")
    q_list <- c(list(p = quantiles), pars)
  }
  thisd <- get_dfun(dist)
  thisq <- get_qfun(dist)



  if(dist == "exp" && pars$par1 <=0) dist <- ""

  thismode <- switch(dist,
                     norm = pars$par1,
                     lnorm = exp(pars$par1 - pars$par2^2),
                     exp = 0,
                     beta = ifelse(test = pars$par1 > 1 & pars$par2 > 1,
                                   yes = (pars$par1-1)/(pars$par1 + pars$par2 - 2),
                                   no = ifelse(test = pars$par1 <=1 & pars$par2 > 1,
                                               yes = 0,
                                               no = ifelse(test = pars$par1 > 1 & pars$par2 <= 1,
                                                           yes = 1,
                                                           no = NA))),
                     gamma = ifelse(pars$par1>1, (pars$par1-1)/pars$par2, NA),
                     weibull = ifelse(test = pars$par1 > 1,
                                      yes = pars$par2*((pars$par1-1)/pars$par1)^(1/pars$par1),
                                      no = 0),
                     do.call(thisq, c(list(optim(par = c(thisval = 0.5),
                                                 lower = c(thisval = 0),
                                                 upper = c(thisval = 1),
                                                 fn = function(thisval, ...) thisd(thisq(thisval, ...), ...) ,
                                                 method = 'L-BFGS-B',
                                                 control = list(fnscale = -1),
                                                 par1 = pars[[1]], par2 = pars[[2]])[[1]][[1]]), pars))
  )

  thismax <- ifelse(!is.na(thismode), do.call(thisd, c(list(thismode), unname(pars))), NA)
  thisquants <- if(is.null(quantiles)) NULL else do.call(thisq, unname(q_list))

  tmp <- c(switch(dist,
                  norm = c(mean = pars$par1,
                           var = pars$par2^2,
                           sd = pars$par2,
                           median = pars$par1),
                  lnorm = c(mean = exp(pars$par1 + 0.5* pars$par2^2),
                            var = (tmpvar <- (exp(2*(pars$par1 + pars$par2^2))-exp(2*pars$par1+pars$par2^2))),
                            sd = sqrt(tmpvar),
                            median = exp(pars$par1)),
                  exp = c(mean = 1/pars$par1,
                          var = (tmpvar <- (1/(pars$par1)^2)),
                          sd = sqrt(tmpvar),
                          median = log(2)/pars$par1),
                  beta = c(mean = (pars$par1/(pars$par1+pars$par2)),
                           var = (tmpvar <- ((pars$par1*pars$par2)/(((pars$par1 + pars$par2)^2)*(pars$par1 + pars$par2 + 1)))),
                           sd = sqrt(tmpvar),
                           median = ifelse(pars$par1 > 1 & pars$par2 > 1, (pars$par1-1/3)/(pars$par1 + pars$par2 -2/3), NA)),
                  gamma = c(mean = pars$par1/pars$par2,
                            var = (tmpvar <- (pars$par1/pars$par2^2)),
                            sd = sqrt(tmpvar),
                            median = NA),
                  weibull = c(mean = pars$par2*gamma(1+1/pars$par1),
                              var = (tmpvar <- ((pars$par2^2)*(gamma(1+2/pars$par1)-(gamma(1+1/pars$par1))^2))),
                              sd = sqrt(tmpvar),
                              median = pars$par2 * (log(2))^(1/pars$par1)),
                  c(mean = NA,
                    var = NA,
                    sd = NA,
                    median = do.call(thisq, c(list(.5), pars)))),
           mode = thismode,
           if(!is.null(quantiles)) thisquants,
           max = thismax)



  if(all(c("mean", "var", "sd", "median", "mode") %in% included)) return(tmp)
  tmp[-which(!c("mean", "var", "sd", "median", "mode") %in% included)]


}
