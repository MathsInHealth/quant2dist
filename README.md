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

When you have distribution parameters for a continuous predictor in different outcome groups, you can estimate effect measures using integration-based methods:

```r
# Example: distribution of biomarker levels by disease status
distdf <- data.frame(
  outcome = c(0, 1),           # 0 = no disease, 1 = disease
  dist = "normal",
  par1 = c(9.3, 8.7),          # means
  par2 = c(1.6, 1.9),          # SDs
  n = c(188, 38)               # sample sizes
)

# Estimate relative risk (RR)
rr_result <- estimate_rr(distdf)
print(rr_result)
#>        RR   log(RR) SE(log(RR))         p    LB2.5   UB97.5
#>  0.931234  -0.07115     0.11234  0.526354 0.747821 1.160234

# Estimate odds ratio (OR)
or_result <- estimate_or(distdf)
print(or_result)
#>        OR   log(OR) SE(log(OR))         p    LB2.5   UB97.5
#>  0.914567  -0.08934     0.13456  0.506789 0.707123 1.183456
```

#### Direct Workflow from Quantiles to OR/RR

Use `q2d_or_rr()` to go directly from quantiles to effect estimates:

```r
# Case group quantiles
case_quants <- c("min" = 5.2, "25%" = 7.1, "50%" = 8.7,
                 "75%" = 10.2, "max" = 15.3, "n" = 38)

# Control group quantiles
control_quants <- c("min" = 6.1, "25%" = 8.2, "50%" = 9.3,
                    "75%" = 10.8, "max" = 14.9, "n" = 188)

# Estimate both OR and RR
result <- q2d_or_rr(case_quants, control_quants, measure = "both")

# View fitted distributions
print(result$distdf)

# View OR result
print(result$OR)

# View RR result
print(result$RR)
```

#### Working with Specific Distribution Families

```r
# Force a specific distribution instead of automatic selection
result_normal <- q2d_or_rr(case_quants, control_quants,
                           measure = "RR",
                           dist_family = "normal")
```

#### Pooling Multiple Studies

`q2d_or_rr()` accepts lists of quantile vectors to pool data from multiple studies:

```r
# Multiple case groups (e.g., from different studies)
cases_list <- list(
  c("min" = 2.1, "25%" = 4.5, "50%" = 6.2, "75%" = 8.3, "max" = 15.2, "n" = 38),
  c("25%" = 3.7, "50%" = 6.0, "75%" = 7.5, "n" = 150)
)

# Multiple control groups
controls_list <- list(
  c("min" = 3.2, "25%" = 5.1, "50%" = 7.1, "75%" = 9.1, "max" = 14.5, "n" = 188),
  c("25%" = 5.5, "50%" = 7.5, "75%" = 9.5, "n" = 200)
)

# Pooled analysis
result_pooled <- q2d_or_rr(cases_list, controls_list)

# distdf will have 4 rows: 2 control groups + 2 case groups
print(result_pooled$distdf)
#>   outcome    dist     par1      par2   n
#> 1       0   lnorm 1.930490 0.3578121 188
#> 2       0 weibull 2.892627 8.4877141 200
#> 3       1   gamma 5.090605 0.7605783  38
#> 4       1    norm 5.719217 2.8151840 150

# Each study can have its own best-fitting distribution
print(result_pooled$RR)
#>        RR   log(RR) SE(log(RR))      p    LB2.5   UB97.5
#>  0.877 -0.132      0.021   <0.001  0.842    0.913
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
cv_risk_rr <- estimate_rr(cv_data)
cv_risk_or <- estimate_or(cv_data)
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

The `estimate_or()` and `estimate_rr()` functions use integration-based maximum likelihood:

**For Odds Ratios** (`estimate_or()`):
1. Fits logistic regression via maximum likelihood
2. Integrates log-likelihood over quantile functions using adaptive quadrature
3. Computes model-based standard errors from the Hessian

**For Relative Risks** (`estimate_rr()`):
1. Uses robust Poisson regression (log-link GLM)
2. Integrates Poisson log-likelihood over quantile functions
3. Applies Zeger-Liang sandwich variance adjustment for robust standard errors

This integration-based approach is more precise than simulation methods, avoiding Monte Carlo error while allowing effect measure estimation when only distributional parameters (not individual-level data) are available.

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
