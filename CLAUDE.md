# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**quant2dist** is an R project for fitting parametric distributions to reported quantile data. It implements maximum likelihood estimation to recover distribution parameters from summary statistics (quantiles, means, SDs) commonly reported in medical/scientific literature. The associated manuscript is "Likelihood estimation for NMA" (network meta-analysis) by Rand et al.

## Refactoring Plan

The goal is to add support for estimation of OR and RR based on estimated quantiles, starting from functions found in the file `R/2026.02.17 - OR_RR_with_example.R`.

These functions should be added to quant2dist, employing similar naming and coding conventions, including CRAN compliance, snake_case for all internal and exported functions, full documentation supporting roxygen2, sensible file splits, etc.



### Step 1: Tidy up and CRAN compliance
- Remove any duplicate code. 
- Identify and fix non-CRAN-compliant patterns: `T`/`F` instead of `TRUE`/`FALSE`, use of `require()` instead of proper imports, `eval(parse(...))` where avoidable, missing function documentation
- Replace `suppressWarnings()` with proper condition handling where appropriate
- Ensure no global variable references or side effects

### Step 2: Rename functions and adopt consistent naming conventions
- Adopt snake_case for all internal and exported functions. 
- Choose clear, descriptive names, rename as seen fit

### Step 3: Organize into a proper R package structure
- Create DESCRIPTION, NAMESPACE (via roxygen2), and standard R package directories
- Keep and ensure maintenance of logically organized R/ source files (e.g. `R/q2d.R`, `R/likelihood.R`, `R/parse.R`, `R/moments.R`, `R/utils.R`)
- Add roxygen2 documentation for all exported functions
- Declare Imports properly 

### Step 4: Test and verify
- Write testthat tests confirming that the new OR/RR function produces identical results to the original
- Test using a mix of distribution families: normal, lognormal, exponential, beta, gamma, weibull
- Run `R CMD check --as-cran` with zero errors, warnings, or notes

### Step 5: Package build and release prep
- Finalize for GitHub: add README.md, .Rbuildignore, LICENSE file
- Ensure `R CMD build` and `R CMD check --as-cran` pass cleanly
- Prepare for eventual CRAN submission

## Current Architecture

The existing quant2dist files are sensibly organized, while the OR/RR estimation function is in a single file entire codebase lives in a single file: `R/2026.02.17 - OR_RR_with_example.R`

### Input format

An example data.frame input is to be found at the bottom of the file.

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
