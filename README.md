# quant2dist

> Fit parametric distributions to reported quantile data and estimate effect measures

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview

**quant2dist** recovers parametric distribution parameters from summary statistics (quantiles, means, standard deviations) commonly reported in medical and scientific literature. The package implements maximum likelihood estimation to fit six distribution families and can estimate odds ratios (OR) and relative risks (RR) from fitted distributions.

### Key Features

- **Distribution fitting** from reported quantiles, means, and SDs
- **Multiple distribution families**: normal, lognormal, exponential, beta, gamma, Weibull
- **Automatic model selection** based on log-likelihood
- **Effect measure estimation**: compute OR and RR from distribution parameters
- **CRAN-compliant** with comprehensive documentation and tests

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("your-username/quant2dist")
```

## Quick Start

### Fitting Distributions to Quantile Data

The main function `q2d()` (quantiles-to-distribution) fits parametric distributions to reported summary statistics:

```r
library(quant2dist)

# Example: reported quantiles from a study
quants <- c(
  "min" = 2.3,
  "25%" = 4.1,
  "50%" = 5.8,  # median
  "75%" = 7.2,
  "max" = 12.1,
  "n" = 50
)

# Fit all distributions and select best
result <- q2d(quants)
result$best_fit
#> $dist
#> [1] "lognormal"
#>
#> $params
#>   meanlog    sdlog
#>   1.6523   0.4201
#>
#> $loglik
#> [1] -87.32
```

### Estimating Odds Ratios and Relative Risks

When you have distribution parameters for a continuous predictor in different outcome groups, you can estimate effect measures:

```r
# Example: distribution of biomarker levels by disease status
# (e.g., from two separate studies or meta-analysis data)
distdf <- data.frame(
  outcome = c(0, 1, 0, 1),           # 0 = no disease, 1 = disease
  dist = "normal",
  par1 = c(9.3, 8.7, 8.4, 7.7),      # means
  par2 = c(1.6, 1.9, 1.3, 1.3),      # SDs
  n = c(188, 38, 188, 38)            # sample sizes
)

# Estimate OR and RR
result <- estimate_or_rr(distdf)
print(result$output)
#>   Type    Ratio logRatio SE(logRatio)       p   LB2.5  UB97.5
#> OR  OR 0.724183 -0.32309     0.143211 0.02413 0.54738 0.95848
#> RR  RR 0.758294 -0.27657     0.126584 0.02891 0.59219 0.97118
```

#### Working with Multiple Distribution Families

```r
# Fit lognormal distributions first
lognormal_data <- data.frame(
  outcome = c(0, 1),
  dist = "lognormal",
  par1 = c(1.5, 1.8),    # meanlog
  par2 = c(0.6, 0.7),    # sdlog
  n = c(100, 50)
)

or_rr_result <- estimate_or_rr(lognormal_data, output = "OR/RR")
print(or_rr_result$output)
```

#### Estimating RR Only

```r
# For cohort studies where RR is preferred
rr_only <- estimate_or_rr(distdf, output = "RR")
print(rr_only$output)
#>   Type    Ratio logRatio SE(logRatio)       p   LB2.5  UB97.5
#> RR  RR 0.758294 -0.27657     0.126584 0.02891 0.59219 0.97118
```

### Using Direct Regression with `rr_glm()`

For binary outcome data with a continuous predictor:

```r
# Example data
dat <- data.frame(
  outcome = c(1, 1, 0, 0, 1, 0, 0, 1),
  biomarker = c(2.1, 3.2, 1.5, 2.0, 4.1, 1.3, 1.8, 3.5)
)

# Fit RR regression with sandwich-adjusted standard errors
rr_coefs <- rr_glm(outcome ~ biomarker, data = dat)
print(rr_coefs)
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept) -0.524221   0.891234 -0.5882   0.5564
#> biomarker    0.312445   0.284512  1.0982   0.2721
```

## Use Cases

### 1. Meta-Analysis and Network Meta-Analysis

Combine studies reporting different summary statistics:

```r
# Study 1 reports: median, IQR
study1 <- c("25%" = 3.2, "50%" = 4.5, "75%" = 6.1, "n" = 45)

# Study 2 reports: mean, SD
study2 <- c("mean" = 4.8, "sd" = 1.9, "n" = 52)

# Fit both and compare
fit1 <- q2d(study1)
fit2 <- q2d(study2)

# Now can combine for OR/RR estimation across studies
```

### 2. Risk Factor Analysis

Estimate effect sizes from published summary statistics:

```r
# Published data: biomarker distribution by cardiovascular event status
cv_data <- data.frame(
  outcome = c(0, 1),  # 0 = no event, 1 = event
  dist = "normal",
  par1 = c(5.2, 6.8),  # mean cholesterol (mmol/L)
  par2 = c(1.1, 1.4),  # SD
  n = c(500, 75)
)

# Estimate per-unit increase in risk
cv_risk <- estimate_or_rr(cv_data)
```

### 3. Simulation and Sample Size Planning

Use fitted distributions for power calculations:

```r
# Fit distribution from pilot data
pilot <- c("min" = 10, "25%" = 15, "50%" = 20, "75%" = 28, "max" = 45, "n" = 30)
pilot_fit <- q2d(pilot)

# Use parameters for simulation
params <- pilot_fit$best_fit$params
# ... simulate data for power analysis
```

## Methodology

### Distribution Fitting

The `q2d()` function uses maximum likelihood estimation to fit six parametric distributions:

- **Normal**: `N(μ, σ²)`
- **Lognormal**: `LN(μ, σ²)`
- **Exponential**: `Exp(λ)`
- **Beta**: `Beta(α, β)`
- **Gamma**: `Gamma(k, θ)`
- **Weibull**: `Weibull(k, λ)`

Model selection is based on maximum log-likelihood across all distributions.

### OR/RR Estimation

The `estimate_or_rr()` function:

1. Generates pseudo-data by sampling from fitted distributions
2. Fits weighted logistic regression for OR estimation
3. Uses case duplication method for RR estimation
4. Applies Zeger-Liang sandwich variance adjustment for robust standard errors

This approach allows estimation of effect measures when only distributional parameters (not individual-level data) are available.

## Citation

If you use quant2dist in your research, please cite:

> Rand D, et al. (2026). Likelihood estimation for network meta-analysis. *Manuscript in preparation*.

## Related Packages

- **estmeansd**: Estimate mean and SD from quantiles (McGrath et al.)
- **metafor**: Meta-analysis framework
- **metamedian**: Meta-analysis of medians

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## Authors

- Daniel Rand (Maths in Health)

---

**Note**: This package is designed for research purposes. Effect size estimates should be interpreted carefully and validated against individual-level data when available.
