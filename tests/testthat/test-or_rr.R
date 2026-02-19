# Tests for integration-based OR/RR estimation functions

# Test data from original implementation
test_distdf <- structure(
  list(n = c(188, 38, 188, 38),
       outcome = c(0, 1, 0, 1),
       dist = c("normal","normal", "normal", "normal"),
       par1 = c(9.3, 8.7, 8.4, 7.7),
       par2 = c(1.6, 1.9, 1.3, 1.3)),
  row.names = 1:4, class = "data.frame")

# Test quantile data for wrapper function
test_case_quants <- c("min" = 5.2, "25%" = 7.1, "50%" = 8.7,
                      "75%" = 10.2, "max" = 15.3, "n" = 38)
test_control_quants <- c("min" = 6.1, "25%" = 8.2, "50%" = 9.3,
                         "75%" = 10.8, "max" = 14.9, "n" = 188)

# =============================================================================
# Tests for estimate_rr()
# =============================================================================

test_that("estimate_rr runs without error on normal distribution", {
  expect_no_error(estimate_rr(test_distdf))
})

test_that("estimate_rr returns correct output structure", {
  result <- estimate_rr(test_distdf)

  # Check it's a named numeric vector
  expect_type(result, "double")
  expect_true(!is.null(names(result)))

  # Check required names
  expect_true(all(c("RR", "log(RR)", "SE(log(RR))", "p", "LB2.5", "UB97.5") %in% names(result)))

  # Check RR is positive
  expect_true(result["RR"] > 0)

  # Check has complete attribute
  expect_true(!is.null(attr(result, "complete")))
  complete <- attr(result, "complete")
  expect_true("vcov" %in% names(complete))
  expect_true("vcov_naive" %in% names(complete))
})

test_that("estimate_rr sandwich SE differs from naive SE", {
  result <- estimate_rr(test_distdf)
  complete <- attr(result, "complete")

  # Sandwich and naive SEs should generally differ
  expect_true("SE_naive" %in% names(complete))
  expect_true(length(complete$SE_naive) == 2)

  # Both should be positive
  expect_true(all(complete$SE > 0))
  expect_true(all(complete$SE_naive > 0))
})

test_that("estimate_rr confidence interval contains point estimate", {
  result <- estimate_rr(test_distdf)

  expect_true(result["RR"] >= result["LB2.5"])
  expect_true(result["RR"] <= result["UB97.5"])
})

test_that("estimate_rr p-value is between 0 and 1", {
  result <- estimate_rr(test_distdf)

  expect_true(result["p"] >= 0)
  expect_true(result["p"] <= 1)
})

test_that("estimate_rr works with lognormal distribution", {
  test_lognormal <- test_distdf
  test_lognormal$dist <- "lognormal"
  test_lognormal$par1 <- log(test_distdf$par1)  # meanlog
  test_lognormal$par2 <- test_distdf$par2 / test_distdf$par1  # approximate sdlog

  result <- expect_no_error(estimate_rr(test_lognormal))
  expect_true(result["RR"] > 0)
})

test_that("estimate_rr works with gamma distribution", {
  test_gamma <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "gamma",
    par1 = c(4, 3.5, 4.2, 3.8),  # shape
    par2 = c(2, 2.2, 1.9, 2.1),  # scale (R parameterization)
    n = c(100, 50, 100, 50)
  )

  result <- expect_no_error(estimate_rr(test_gamma))
  expect_true(result["RR"] > 0)
})

# =============================================================================
# Tests for estimate_or()
# =============================================================================

test_that("estimate_or runs without error on normal distribution", {
  expect_no_error(estimate_or(test_distdf))
})

test_that("estimate_or returns correct output structure", {
  result <- estimate_or(test_distdf)

  # Check it's a named numeric vector
  expect_type(result, "double")
  expect_true(!is.null(names(result)))

  # Check required names
  expect_true(all(c("OR", "log(OR)", "SE(log(OR))", "p", "LB2.5", "UB97.5") %in% names(result)))

  # Check OR is positive
  expect_true(result["OR"] > 0)

  # Check has complete attribute
  expect_true(!is.null(attr(result, "complete")))
})

test_that("estimate_or confidence interval contains point estimate", {
  result <- estimate_or(test_distdf)

  expect_true(result["OR"] >= result["LB2.5"])
  expect_true(result["OR"] <= result["UB97.5"])
})

test_that("estimate_or p-value is between 0 and 1", {
  result <- estimate_or(test_distdf)

  expect_true(result["p"] >= 0)
  expect_true(result["p"] <= 1)
})

test_that("estimate_or works with weibull distribution", {
  test_weibull <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "weibull",
    par1 = c(2.5, 2.3, 2.6, 2.4),  # shape
    par2 = c(8, 7, 8.5, 7.5),      # scale
    n = c(100, 50, 100, 50)
  )

  result <- expect_no_error(estimate_or(test_weibull))
  expect_true(result["OR"] > 0)
})

test_that("estimate_or works with beta distribution", {
  # Beta requires values in (0,1), so scale parameters accordingly
  test_beta <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "beta",
    par1 = c(3, 2.5, 3.2, 2.8),  # shape1
    par2 = c(5, 4.5, 5.2, 4.8),  # shape2
    n = c(100, 50, 100, 50)
  )

  result <- expect_no_error(estimate_or(test_beta))
  expect_true(result["OR"] > 0)
})

# =============================================================================
# Tests for q2d_or_rr() wrapper function
# =============================================================================

test_that("q2d_or_rr runs without error", {
  expect_no_error(q2d_or_rr(test_case_quants, test_control_quants))
})

test_that("q2d_or_rr returns correct structure", {
  result <- q2d_or_rr(test_case_quants, test_control_quants)

  # Check list components
  expect_type(result, "list")
  expect_true(all(c("case_fit", "control_fit", "distdf", "OR", "RR") %in% names(result)))

  # Check distdf structure
  expect_s3_class(result$distdf, "data.frame")
  expect_equal(nrow(result$distdf), 2)
  expect_equal(result$distdf$outcome, c(0, 1))

  # Check OR and RR present
  expect_true(!is.null(result$OR))
  expect_true(!is.null(result$RR))
})

test_that("q2d_or_rr respects measure parameter", {
  # OR only
  result_or <- q2d_or_rr(test_case_quants, test_control_quants, measure = "OR")
  expect_true("OR" %in% names(result_or))
  expect_false("RR" %in% names(result_or))

  # RR only
  result_rr <- q2d_or_rr(test_case_quants, test_control_quants, measure = "RR")
  expect_false("OR" %in% names(result_rr))
  expect_true("RR" %in% names(result_rr))
})

test_that("q2d_or_rr respects dist_family parameter", {
  result <- q2d_or_rr(test_case_quants, test_control_quants, dist_family = "normal")

  # Both should be normal
  expect_equal(result$distdf$dist[1], "norm")
  expect_equal(result$distdf$dist[2], "norm")
})

test_that("q2d_or_rr errors on missing n", {
  bad_quants <- test_case_quants[names(test_case_quants) != "n"]

  expect_error(
    q2d_or_rr(bad_quants, test_control_quants),
    "must include 'n'"
  )
})

test_that("q2d_or_rr errors on invalid measure", {
  expect_error(
    q2d_or_rr(test_case_quants, test_control_quants, measure = "invalid"),
    "must be 'OR', 'RR', or 'both'"
  )
})

# =============================================================================
# Tests for list inputs (multiple studies/groups)
# =============================================================================

test_that("q2d_or_rr accepts list of case vectors", {
  cases_list <- list(
    c("min" = 5.2, "25%" = 7.1, "50%" = 8.7, "75%" = 10.2, "max" = 15.3, "n" = 38),
    c("25%" = 6.5, "50%" = 8.0, "75%" = 9.5, "n" = 50)
  )

  result <- expect_no_error(q2d_or_rr(cases_list, test_control_quants, measure = "RR"))

  # Should have 3 rows: 1 control + 2 case groups
  expect_equal(nrow(result$distdf), 3)
  expect_equal(sum(result$distdf$outcome == 0), 1)
  expect_equal(sum(result$distdf$outcome == 1), 2)

  # case_fit should be a list
  expect_true(is.list(result$case_fit))
  expect_equal(length(result$case_fit), 2)
})

test_that("q2d_or_rr accepts list of control vectors", {
  controls_list <- list(
    c("min" = 6.1, "25%" = 8.2, "50%" = 9.3, "75%" = 10.8, "max" = 14.9, "n" = 188),
    c("25%" = 8.0, "50%" = 9.5, "75%" = 11.0, "n" = 200)
  )

  result <- expect_no_error(q2d_or_rr(test_case_quants, controls_list, measure = "OR"))

  # Should have 3 rows: 2 control + 1 case group
  expect_equal(nrow(result$distdf), 3)
  expect_equal(sum(result$distdf$outcome == 0), 2)
  expect_equal(sum(result$distdf$outcome == 1), 1)

  # control_fit should be a list
  expect_true(is.list(result$control_fit))
  expect_equal(length(result$control_fit), 2)
})

test_that("q2d_or_rr accepts lists for both cases and controls", {
  cases_list <- list(
    c("25%" = 6.5, "50%" = 8.0, "75%" = 9.5, "n" = 50),
    c("25%" = 7.0, "50%" = 8.5, "75%" = 10.0, "n" = 60)
  )

  controls_list <- list(
    c("25%" = 8.0, "50%" = 9.5, "75%" = 11.0, "n" = 200),
    c("25%" = 8.5, "50%" = 10.0, "75%" = 11.5, "n" = 250)
  )

  result <- expect_no_error(q2d_or_rr(cases_list, controls_list, measure = "both"))

  # Should have 4 rows: 2 control + 2 case groups
  expect_equal(nrow(result$distdf), 4)
  expect_equal(sum(result$distdf$outcome == 0), 2)
  expect_equal(sum(result$distdf$outcome == 1), 2)

  # Both should be lists
  expect_true(is.list(result$case_fit))
  expect_true(is.list(result$control_fit))

  # Should have both OR and RR
  expect_true("OR" %in% names(result))
  expect_true("RR" %in% names(result))
})

test_that("q2d_or_rr list inputs produce valid effect estimates", {
  cases_list <- list(
    c("25%" = 6.5, "50%" = 8.0, "75%" = 9.5, "n" = 50),
    c("25%" = 7.0, "50%" = 8.5, "75%" = 10.0, "n" = 60)
  )

  controls_list <- list(
    c("25%" = 8.0, "50%" = 9.5, "75%" = 11.0, "n" = 200),
    c("25%" = 8.5, "50%" = 10.0, "75%" = 11.5, "n" = 250)
  )

  result <- q2d_or_rr(cases_list, controls_list)

  # Check OR is positive
  expect_true(result$OR["OR"] > 0)
  expect_true(result$OR["SE(log(OR))"] > 0)
  expect_true(result$OR["p"] >= 0 && result$OR["p"] <= 1)

  # Check RR is positive
  expect_true(result$RR["RR"] > 0)
  expect_true(result$RR["SE(log(RR))"] > 0)
  expect_true(result$RR["p"] >= 0 && result$RR["p"] <= 1)
})

test_that("q2d_or_rr errors on list with missing n", {
  bad_list <- list(
    c("25%" = 6.5, "50%" = 8.0, "75%" = 9.5, "n" = 50),
    c("25%" = 7.0, "50%" = 8.5, "75%" = 10.0)  # missing n
  )

  expect_error(
    q2d_or_rr(bad_list, test_control_quants),
    "must include 'n'"
  )
})

# =============================================================================
# Tests for consistency and edge cases
# =============================================================================

test_that("OR and RR are similar for rare outcomes", {
  # For rare outcomes, OR should approximate RR
  rare_outcome <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "normal",
    par1 = c(10, 9.5, 10.2, 9.7),
    par2 = c(1.5, 1.6, 1.4, 1.5),
    n = c(1000, 10, 1000, 10)  # ~1% prevalence
  )

  or_result <- estimate_or(rare_outcome)
  rr_result <- estimate_rr(rare_outcome)

  # OR and RR should be similar (within 10%)
  expect_true(abs(or_result["OR"] - rr_result["RR"]) / rr_result["RR"] < 0.15)
})

test_that("Results are consistent with same data", {
  # Running twice should give identical results
  result1 <- estimate_rr(test_distdf)
  result2 <- estimate_rr(test_distdf)

  expect_equal(result1, result2)
})

test_that("estimate_or handles custom column names", {
  custom_df <- test_distdf
  colnames(custom_df) <- c("sample_size", "event", "distribution", "param1", "param2")

  result <- estimate_or(
    custom_df,
    colnms = c(n = "sample_size", outcome = "event", dist = "distribution",
               par1 = "param1", par2 = "param2")
  )

  expect_true(result["OR"] > 0)
})

test_that("estimate_rr handles custom column names", {
  custom_df <- test_distdf
  colnames(custom_df) <- c("sample_size", "event", "distribution", "param1", "param2")

  result <- estimate_rr(
    custom_df,
    colnms = c(n = "sample_size", outcome = "event", dist = "distribution",
               par1 = "param1", par2 = "param2")
  )

  expect_true(result["RR"] > 0)
})

test_that("Functions error appropriately on missing columns", {
  bad_df <- test_distdf[, -1]  # Remove n column

  expect_error(estimate_or(bad_df), "does not contain necessary information")
  expect_error(estimate_rr(bad_df), "does not contain necessary information")
})

# =============================================================================
# Tests for numerical stability
# =============================================================================

test_that("Functions handle extreme parameter values", {
  extreme_df <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "normal",
    par1 = c(100, 95, 102, 97),   # Large values
    par2 = c(0.1, 0.15, 0.08, 0.12),  # Small SDs
    n = c(50, 50, 50, 50)
  )

  expect_no_error(estimate_or(extreme_df))
  expect_no_error(estimate_rr(extreme_df))
})

test_that("Standard errors are always positive", {
  result_or <- estimate_or(test_distdf)
  result_rr <- estimate_rr(test_distdf)

  expect_true(result_or["SE(log(OR))"] > 0)
  expect_true(result_rr["SE(log(RR))"] > 0)

  # Also check complete output
  expect_true(all(attr(result_rr, "complete")$SE > 0))
  expect_true(all(attr(result_rr, "complete")$SE_naive > 0))
})
