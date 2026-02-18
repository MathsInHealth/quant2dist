#' Compute robust sandwich variance estimator
#'
#' Implements the Zeger-Liang sandwich variance estimator for GLM models
#' to account for within-cluster correlation. Used internally for
#' adjusting standard errors in relative risk estimation.
#'
#' @param m A fitted GLM object (from \code{glm()}).
#' @param unique_val Character string naming the clustering variable in
#'   \code{m$data}. Default is \code{"id"}.
#'
#' @return A list with four components: \code{bread} (model-based variance-covariance matrix),
#'   \code{meat} (empirical cross-product sum), \code{sandwich} (robust variance-covariance matrix),
#'   and \code{SE} (vector of robust standard errors).
#'
#' @references
#'   Zeger SL, Liang KY. Longitudinal data analysis for discrete and continuous outcomes.
#'   Biometrics. 1986;42(1):121-130.
#'
#' @importFrom stats residuals model.matrix
#' @keywords internal
zl_sandwich <- function(m, unique_val = 'id') {
  mres <- residuals(m, 'response')*sqrt(m$prior.weights)
  mmat <- model.matrix(m)
  sm <- summary(m)
  bread <- sm$cov.unscaled
  meal <- lapply(unique(m$data[,unique_val]), function(i) {
    irows <- which(m$data[,unique_val] == i)
    crossprod(mres[irows] %*% mmat[irows,,drop=F])
  })
  meat <- Reduce('+', meal)
  sandwich <- bread %*% meat %*% bread
  list(bread = bread, meat = meat, sandwich = sandwich, SE = sqrt(diag(sandwich)))
}

#' Estimate odds ratios and relative risks from distribution parameters
#'
#' Computes odds ratios (OR) and relative risks (RR) relating a continuous
#' predictor to a binary outcome, given parametric distributions fit to the
#' predictor in each outcome group. Uses simulation-based pseudo-data with
#' weighted logistic regression. Standard errors for RR are adjusted using
#' the Zeger-Liang sandwich estimator.
#'
#' @param distdf A data frame with one row per outcome group, containing:
#'   \itemize{
#'     \item \code{outcome}: Binary outcome indicator (0 or 1).
#'     \item \code{dist}: Distribution name (e.g., \code{"normal"}, \code{"lognormal"}).
#'     \item \code{par1}: First distribution parameter (e.g., mean or shape).
#'     \item \code{par2}: Second distribution parameter (e.g., SD or scale).
#'     \item \code{n}: Sample size for each group.
#'   }
#' @param colnms Named character vector mapping standard column names to
#'   actual column names in \code{distdf}. Default assumes standard names.
#' @param stvs Starting values (currently unused; reserved for future extensions).
#' @param nsamples Number of pseudo-observations to simulate per distribution.
#'   Default is 1000.
#' @param output Character string controlling which estimates to compute:
#'   \code{"OR"} for odds ratio only, \code{"RR"} for relative risk only,
#'   or \code{"OR/RR"} (default) for both.
#'
#' @return A list containing: \code{row_output} (named numeric vector of estimates),
#'   \code{output} (data frame with one row per measure containing ratio, log-ratio, SE, p-value, and CIs),
#'   and when computed: \code{OR_coef}, \code{OR_coefficients}, \code{OR_glm} (OR model components),
#'   \code{RR_coef}, \code{RR_coefficients}, \code{RR_glm}, \code{RR_vcov} (RR model components with
#'   sandwich-adjusted standard errors).
#'
#' @details
#' For OR estimation, standard logistic regression is used on simulated data.
#' For RR estimation, the method duplicates cases with \code{outcome = 1},
#' setting the duplicates to \code{outcome = 0}, then applies logistic regression
#' with sandwich variance adjustment to approximate relative risk.
#'
#' Supported distributions include: normal, lognormal, exponential, beta, gamma, weibull.
#'
#' @examples
#' \dontrun{
#' # Example: estimate OR/RR for normal distributions
#' df <- data.frame(
#'   outcome = c(0, 1, 0, 1),
#'   dist = "normal",
#'   par1 = c(9.3, 8.7, 8.4, 7.7),
#'   par2 = c(1.6, 1.9, 1.3, 1.3),
#'   n = c(188, 38, 188, 38)
#' )
#' result <- estimate_or_rr(df)
#' print(result$output)
#' }
#'
#' @importFrom stats glm binomial qnorm pt
#' @export
estimate_or_rr <- function(distdf,
                      colnms = c(outcome = 'outcome', dist = 'dist', par1 = 'par1', par2 = 'par2', n = 'n'),
                      stvs = NULL,
                      nsamples = 1000,
                      output = 'OR/RR') {

  do_OR <- length(grep(pattern = 'OR', fixed = TRUE, x = toupper(output)))
  do_RR <- length(grep(pattern = 'RR', fixed = TRUE, x = toupper(output)))

  # Ensuring that correct coumn information is available and correctly named
  thesecolnms <- c(outcome = 'outcome', dist = 'dist', par1 = 'par1', par2 = 'par2', n = 'n')
  colnms <- colnms[names(colnms) %in% names(thesecolnms)]
  if(length(colnms)) thesecolnms[names(colnms)] <- colnms

  if(!all(thesecolnms %in% colnames(distdf))) stop('distdf does not contain necessary information.')

  distdf <- distdf[, thesecolnms]
  colnames(distdf) <- names(thesecolnms)

  # Generate simulated data for OR estimation
  sim_df <- do.call(rbind, lapply(1:NROW(distdf), function(i) {
    data.frame(outcome = distdf$outcome[i],
               measure = do.call(get_qfun(distdf$dist[i]), list((1:nsamples-.5)/nsamples, distdf$par1[i], distdf$par2[i])),
               weight = distdf$n[i]/nsamples)
  }))
  sim_df$id <- 1:NROW(sim_df)

  # OR estimation
  # Note: suppressWarnings used to handle GLM convergence warnings from weighted pseudodata
  or_glm <- suppressWarnings(glm(outcome ~ measure,
                                 data = sim_df,
                                 weights = sim_df$weight,
                                 family = binomial(link = 'logit')))
  if(do_RR) {
    # Generate simulated data for RR estimation
    sim_df2 <- sim_df[which(sim_df$outcome ==1),]
    sim_df2$outcome <- 0
    sim_df2 <- rbind(sim_df, sim_df2)
    sim_df2 <- sim_df2[order(sim_df2$id),]

    # RR estimation
    # Note: suppressWarnings used to handle GLM convergence warnings from weighted pseudodata
    rr_glm <- suppressWarnings(glm(outcome ~ measure,
                                   data = sim_df2,
                                   weights = sim_df2$weight,
                                   family = binomial(link = 'logit')))
    # summary(rr_glm)$coefficients

    # Standard error correction using Zeber Liang sandwich
    ZL <- zl_sandwich(rr_glm)
    rrc <- summary(rr_glm)$coefficients
    # Corrected SE
    rrc[,2] <- unname(ZL$SE)
    rrc[,3] <- abs(rrc[,1])/rrc[,2]
    rrc[,4] <- 2*pt(rrc[,3], sum(sim_df2$weight)-1, lower.tail = FALSE)
  }


  orc <- summary(or_glm)$coefficients

  # Note: uses stats::qnorm - ensure proper import
  rout <- c(
    if(do_OR)
      c(OR =            exp(orc[2,1]),
        'log(OR)' =     orc[2,1],
        'SE(log(OR))' = orc[2,2],
        'p(OR)' =       orc[2,4],
        'LB2.5(OR)' =   exp(orc[2,1] + qnorm(0.025)*orc[2,2]),
        'UB97.5(OR)' =  exp(orc[2,1] + qnorm(0.975)*orc[2,2])) else NULL,
    if(do_RR)
      c(RR =            exp(rrc[2,1]),
        'log(RR)' =     rrc[2,1],
        'SE(log(RR))' = rrc[2,2],
        'p(RR)' =       rrc[2,4],
        'LB2.5(RR)' =   exp(rrc[2,1] + qnorm(0.025)*rrc[2,2]),
        'UB97.5(RR)' =  exp(rrc[2,1] + qnorm(0.975)*rrc[2,2])) else NULL
  )
  outp <- as.data.frame(rbind(
    if(do_OR) rout[grep(pattern = 'OR', x = names(rout), value = FALSE, fixed = TRUE)] else NULL,
    if(do_RR) rout[grep(pattern = 'RR', x = names(rout), value = FALSE, fixed = TRUE)] else NULL
  ))
  outp <- cbind(c(if(do_OR) 'OR' else NULL,
                  if(do_RR) 'RR' else NULL),
                outp)

  colnames(outp) <- c('Type', 'Ratio', 'logRatio', 'SE(logRatio)', 'p', 'LB2.5', 'UB97.5')

  rownames(outp) <- outp$Type
  c(list(
    row_output = rout,
    output = outp),
    if(do_OR) list(
      OR_coef = or_glm$coefficients,
      OR_coefficients = summary(or_glm)$coefficients,
      OR_glm = or_glm) else NULL,
    if(do_RR) list(
      RR_coef = rr_glm$coefficients,
      RR_coefficients = rrc,
      RR_glm = rr_glm,
      RR_vcov = ZL$sandwich) else NULL)


}

#' Fit a relative risk regression model using logistic approximation
#'
#' Estimates relative risks (rather than odds ratios) from binary outcome data
#' using the case duplication method with sandwich variance adjustment.
#' Standard logistic regression on duplicated data approximates the RR,
#' and robust standard errors account for the induced correlation.
#'
#' @param formula A formula object specifying the model (e.g., \code{outcome ~ predictor}).
#' @param data A data frame containing the variables in \code{formula}.
#' @param weights Numeric vector of case weights, or a single value (default 1)
#'   applied to all observations.
#'
#' @return A matrix containing the coefficient table with sandwich-adjusted
#'   standard errors, z-statistics, and p-values. Rows correspond to model
#'   terms (intercept, predictors); columns are Estimate, Std. Error, z value, Pr(>|z|).
#'
#' @details
#' This function duplicates all cases where the outcome equals 1, setting
#' the outcome to 0 in the duplicates, then fits a logistic regression model.
#' The resulting log-odds coefficient approximates log(RR). Standard errors
#' are corrected using \code{\link{zl_sandwich}} to account for the case duplication.
#'
#' An \code{id} column is added to \code{data} if not already present to track
#' original vs. duplicated observations.
#'
#' @examples
#' \dontrun{
#' # Simple RR regression
#' dat <- data.frame(outcome = c(1,1,0,0,1,0), x = c(2,3,1,2,4,1))
#' rr_glm(outcome ~ x, data = dat)
#' }
#'
#' @importFrom stats glm binomial pt
#' @export
rr_glm <- function(formula, data, weights = 1) {
  # Generate simulated data for RR estimation
  outv <- as.character(formula[[2]])
  # data$outcome <- data[,outv]

  # Create copies to avoid modifying input data
  if(!'id' %in% colnames(data)) {
    data$id <- 1:NROW(data)
  }

  data$weights <- weights

  data2 <- data[which(data[, outv] ==1),]
  data2[, outv] <- 0
  data2 <- rbind(data, data2)
  data2 <- data2[order(data2$id),]

  # RR estimation
  # Note: suppressWarnings used to handle GLM convergence warnings from case duplication method
  rr_glm <- suppressWarnings(do.call(glm, list(formula = formula,
                                 data = data2,
                                 weights = data2$weights,
                                 family = binomial(link = 'logit'))))
  # summary(rr_glm)$coefficients

  # Standard error correction using Zeber Liang sandwich
  ZL <- zl_sandwich(rr_glm)
  rrc <- summary(rr_glm)$coefficients
  # Corrected SE
  rrc[,2] <- unname(ZL$SE)
  rrc[,3] <- abs(rrc[,1])/rrc[,2]
  rrc[,4] <- 2*pt(rrc[,3], sum(data2$weights)-1, lower.tail = FALSE)
  rrc
}
