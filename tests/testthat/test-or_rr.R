# Tests for OR/RR estimation functions

# Test data from original implementation
test_input <- structure(
  list(n = c(188, 38, 188, 38),
       outcome = c(0, 1, 0, 1),
       dist = c("normal","normal", "normal", "normal"),
       par1 = c(9.3, 8.7, 8.4, 7.7),
       par2 = c(1.6, 1.9, 1.3, 1.3)),
  row.names = 1:4, class = "data.frame")

test_that("estimate_or_rr runs without error on normal distribution", {
  expect_no_error(estimate_or_rr(test_input))
})

test_that("estimate_or_rr returns correct output structure", {
  result <- estimate_or_rr(test_input)

  # Check list components exist
  expect_type(result, "list")
  expect_true("row_output" %in% names(result))
  expect_true("output" %in% names(result))
  expect_true("OR_coef" %in% names(result))
  expect_true("OR_coefficients" %in% names(result))
  expect_true("OR_glm" %in% names(result))
  expect_true("RR_coef" %in% names(result))
  expect_true("RR_coefficients" %in% names(result))
  expect_true("RR_glm" %in% names(result))
  expect_true("RR_vcov" %in% names(result))

  # Check output data frame structure
  expect_s3_class(result$output, "data.frame")
  expect_equal(nrow(result$output), 2)  # OR and RR rows
  expect_equal(result$output$Type, c("OR", "RR"))
  expect_true(all(c("Type", "Ratio", "logRatio", "SE(logRatio)", "p", "LB2.5", "UB97.5") %in% colnames(result$output)))

  # Check that estimates are numeric and positive
  expect_true(is.numeric(result$output$Ratio))
  expect_true(all(result$output$Ratio > 0))
})

test_that("estimate_or_rr returns OR only when requested", {
  result <- estimate_or_rr(test_input, output = "OR")

  # Should have OR components
  expect_true("OR_coef" %in% names(result))
  expect_true("OR_coefficients" %in% names(result))
  expect_true("OR_glm" %in% names(result))

  # Should not have RR components
  expect_false("RR_coef" %in% names(result))
  expect_false("RR_coefficients" %in% names(result))
  expect_false("RR_glm" %in% names(result))
  expect_false("RR_vcov" %in% names(result))

  # Output should have only OR row
  expect_equal(nrow(result$output), 1)
  expect_equal(result$output$Type, "OR")
})

test_that("estimate_or_rr returns RR only when requested", {
  result <- estimate_or_rr(test_input, output = "RR")

  # Should have RR components
  expect_true("RR_coef" %in% names(result))
  expect_true("RR_coefficients" %in% names(result))
  expect_true("RR_glm" %in% names(result))
  expect_true("RR_vcov" %in% names(result))

  # Should not have OR components
  expect_false("OR_coef" %in% names(result))
  expect_false("OR_coefficients" %in% names(result))
  expect_false("OR_glm" %in% names(result))

  # Output should have only RR row
  expect_equal(nrow(result$output), 1)
  expect_equal(result$output$Type, "RR")
})

test_that("estimate_or_rr works with lognormal distribution", {
  test_lognormal <- test_input
  test_lognormal$dist <- "lognormal"

  result <- expect_no_error(estimate_or_rr(test_lognormal))
  expect_s3_class(result$output, "data.frame")
  expect_equal(nrow(result$output), 2)
})

test_that("estimate_or_rr works with exponential distribution", {
  test_exp <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "exponential",
    par1 = c(0.5, 0.3, 0.6, 0.4),  # rate parameter
    par2 = c(1, 1, 1, 1),  # exponential only uses par1, but need par2
    n = c(100, 50, 100, 50)
  )

  result <- expect_no_error(estimate_or_rr(test_exp))
  expect_s3_class(result$output, "data.frame")
})

test_that("estimate_or_rr works with gamma distribution", {
  test_gamma <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "gamma",
    par1 = c(2, 2.5, 1.8, 2.2),  # shape
    par2 = c(0.5, 0.6, 0.4, 0.55),  # rate
    n = c(100, 50, 100, 50)
  )

  result <- expect_no_error(estimate_or_rr(test_gamma))
  expect_s3_class(result$output, "data.frame")
})

test_that("estimate_or_rr works with weibull distribution", {
  test_weibull <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "weibull",
    par1 = c(2, 2.2, 1.9, 2.1),  # shape
    par2 = c(5, 4.5, 5.2, 4.8),  # scale
    n = c(100, 50, 100, 50)
  )

  result <- expect_no_error(estimate_or_rr(test_weibull))
  expect_s3_class(result$output, "data.frame")
})

test_that("estimate_or_rr works with beta distribution", {
  test_beta <- data.frame(
    outcome = c(0, 1, 0, 1),
    dist = "beta",
    par1 = c(2, 2.5, 1.8, 2.3),  # shape1
    par2 = c(3, 3.5, 2.8, 3.2),  # shape2
    n = c(100, 50, 100, 50)
  )

  result <- expect_no_error(estimate_or_rr(test_beta))
  expect_s3_class(result$output, "data.frame")
})

test_that("estimate_or_rr respects nsamples parameter", {
  # Smaller nsamples should run faster but be less precise
  result_small <- estimate_or_rr(test_input, nsamples = 100)
  result_large <- estimate_or_rr(test_input, nsamples = 2000)

  # Both should return valid results
  expect_s3_class(result_small$output, "data.frame")
  expect_s3_class(result_large$output, "data.frame")

  # Results should be similar but not identical
  expect_true(abs(result_small$output$Ratio[1] - result_large$output$Ratio[1]) < 0.5)
})

test_that("estimate_or_rr handles custom column names", {
  test_custom <- test_input
  colnames(test_custom) <- c("sample_size", "event", "distribution", "param1", "param2")

  result <- estimate_or_rr(
    test_custom,
    colnms = c(n = "sample_size", outcome = "event", dist = "distribution",
               par1 = "param1", par2 = "param2")
  )

  expect_s3_class(result$output, "data.frame")
  expect_equal(nrow(result$output), 2)
})

test_that("estimate_or_rr errors on missing columns", {
  test_missing <- test_input[, -1]  # Remove n column

  expect_error(
    estimate_or_rr(test_missing),
    "does not contain necessary information"
  )
})

test_that("rr_glm runs without error", {
  dat <- data.frame(
    outcome = c(1, 1, 0, 0, 1, 0, 0, 1),
    x = c(2, 3, 1, 2, 4, 1, 1.5, 3.5)
  )

  result <- expect_no_error(rr_glm(outcome ~ x, data = dat))

  # Check it returns a matrix
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4)  # Estimate, Std. Error, z value, Pr(>|z|)
  expect_equal(nrow(result), 2)  # Intercept + x
})

test_that("rr_glm respects weights parameter", {
  dat <- data.frame(
    outcome = c(1, 1, 0, 0, 1, 0),
    x = c(2, 3, 1, 2, 4, 1)
  )

  result_unweighted <- rr_glm(outcome ~ x, data = dat, weights = 1)
  result_weighted <- rr_glm(outcome ~ x, data = dat, weights = c(1, 1, 1, 1, 2, 1))

  # Both should run successfully
  expect_true(is.matrix(result_unweighted))
  expect_true(is.matrix(result_weighted))

  # Results should differ
  expect_false(identical(result_unweighted, result_weighted))
})

test_that("zl_sandwich returns correct structure", {
  # Create a simple GLM to test sandwich estimator
  dat <- data.frame(
    y = c(1, 0, 1, 0, 1, 1, 0, 0),
    x = c(1, 2, 3, 4, 5, 6, 7, 8),
    id = c(1, 1, 2, 2, 3, 3, 4, 4)
  )

  m <- glm(y ~ x, data = dat, family = binomial())

  result <- quant2dist:::zl_sandwich(m, unique_val = 'id')

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("bread", "meat", "sandwich", "SE") %in% names(result)))

  # Check dimensions
  expect_true(is.matrix(result$bread))
  expect_true(is.matrix(result$meat))
  expect_true(is.matrix(result$sandwich))
  expect_true(is.numeric(result$SE))
  expect_equal(length(result$SE), 2)  # Intercept + x
})

test_that("Confidence intervals contain point estimates", {
  result <- estimate_or_rr(test_input)

  # For OR
  or_row <- result$output[result$output$Type == "OR", ]
  expect_true(or_row$Ratio >= or_row$LB2.5)
  expect_true(or_row$Ratio <= or_row$UB97.5)

  # For RR
  rr_row <- result$output[result$output$Type == "RR", ]
  expect_true(rr_row$Ratio >= rr_row$LB2.5)
  expect_true(rr_row$Ratio <= rr_row$UB97.5)
})

test_that("Standard errors are positive", {
  result <- estimate_or_rr(test_input)

  expect_true(all(result$output$`SE(logRatio)` > 0))
})

test_that("P-values are between 0 and 1", {
  result <- estimate_or_rr(test_input)

  expect_true(all(result$output$p >= 0))
  expect_true(all(result$output$p <= 1))
})
