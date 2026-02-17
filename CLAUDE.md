# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**quant2dist** is an R project for fitting parametric distributions to reported quantile data. It implements maximum likelihood estimation to recover distribution parameters from summary statistics (quantiles, means, SDs) commonly reported in medical/scientific literature. The associated manuscript is "Likelihood estimation for NMA" (network meta-analysis) by Rand et al.

## Refactoring Plan

The goal is to convert this script-based project into a CRAN-compliant R package named **quant2dist**, with `q2d()` as the primary user-facing function (replacing `qest()`). This should be done in the following sequential steps:

### Step 1: Tidy the source code
- Remove all graph/plot-related functions (`plot_qest`, `table_qest`, `ir_graph`, and any plotting helpers)
- Remove all commented-out code blocks
- Remove obsolete or unused functions (e.g. `contvar2` if superseded by `contvar2s`/`contvar2n`, duplicate helpers)
- Keep only the code paths required by `qest()` and its dependencies

### Step 2: CRAN compliance fixes
- Identify and fix non-CRAN-compliant patterns: `T`/`F` instead of `TRUE`/`FALSE`, use of `require()` instead of proper imports, `eval(parse(...))` where avoidable, missing function documentation
- Replace `suppressWarnings()` with proper condition handling where appropriate
- Ensure no global variable references or side effects

### Step 3: Rename functions and adopt consistent naming conventions
- Rename `qest()` → `q2d()` as the primary entry point
- Adopt snake_case for all internal and exported functions
- Choose clear, descriptive names (e.g. `parse_quantiles()` instead of `idf2()`)

### Step 4: Organize into a proper R package structure
- Create DESCRIPTION, NAMESPACE (via roxygen2), and standard R package directories
- Split the single monolithic file into logically organized R/ source files (e.g. `R/q2d.R`, `R/likelihood.R`, `R/parse.R`, `R/moments.R`, `R/utils.R`)
- Add roxygen2 documentation for all exported functions
- Declare Imports properly (statmod for Gauss-Hermite quadrature; no ggplot2/gt/ggsci dependencies needed after graph removal)

### Step 5: Test and verify
- Write testthat tests confirming `q2d()` produces identical results to the original `qest()`
- Test across all 6 distribution families: normal, lognormal, exponential, beta, gamma, weibull
- Verify edge cases (single quantile, min/max inputs, mean+SD inputs)
- Run `R CMD check --as-cran` with zero errors, warnings, or notes

### Step 6: Package build and release prep
- Finalize for GitHub: add README.md, .Rbuildignore, LICENSE file
- Ensure `R CMD build` and `R CMD check --as-cran` pass cleanly
- Prepare for eventual CRAN submission

## Current Architecture

The entire codebase lives in a single file: `R/2022.09.29 - quantile_fit.R` (~2700 lines). There is no package structure yet.

### Key function layers (pre-refactor)

1. **Likelihood functions** — compute negative log-likelihood for distribution fitting:
   - `contvar2()` — original version, distribution-specific branching
   - `contvar2s()` — streamlined version using generic `dfun`/`pfun` dispatch
   - `contvar2n()` — further refined variant (default used by `qest`)
   - `contvar2r()` — variant with random effects
   - `ir2()`, `ir2r()` — individual-record likelihood functions
   - All accept `dist` parameter: `"normal"`, `"lognormal"`, `"exponential"`, `"beta"`, `"gamma"`, `"weibull"`

2. **Fitting functions**:
   - `ir_fit2()`, `ir_fit3()` — low-level optimizers wrapping `optim()`
   - `qest()` — **main entry point**: fits all 6 distributions, selects best by log-likelihood

3. **Data preparation**:
   - `idf()`, `idf2()` — convert a named vector of quantiles into internal data frame format (`lb`, `ub`, `dtype`, `q`, `r`, `N`)
   - `qefit_sample()` — generates starting values for all distributions via method-of-moments
   - `quants_df()` — converts a data frame of quantile columns into the internal format

4. **Utilities**:
   - `meanSD()`, `meanSE()` — compute moments for each distribution family
   - `dodist()` — generic wrapper for d/p/q/r functions across distributions
   - `fixdist()` — normalize distribution name strings
   - `joint_p()`, `joint_pnorm()` — joint probability calculations for mean/SD observations
   - `moments()`, `moments_df()` — compute distribution moments and quantiles from parameters

### Data type conventions

The `dtype` field in internal data frames controls how observations are treated:
- `dtype == 1` — quantile observation (point or interval)
- `dtype == 2` — observed mean (no SD)
- `dtype == 3` — observed mean AND SD (joint probability via normal product integral)

### Input format

The standard input is a named numeric vector where names are quantile labels (e.g., `"25%"`, `"50%"`, `"min"`, `"max"`) plus `"n"` for sample size. The `idf()`/`idf2()` functions parse these into the internal data structure.

## Build and Check Commands

```bash
# Build the package (once package structure exists)
R CMD build .

# Run CRAN checks
R CMD check --as-cran quant2dist_*.tar.gz

# Run tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Generate documentation from roxygen2
Rscript -e 'roxygen2::roxygenise()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-q2d.R")'
```
