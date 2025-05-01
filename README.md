<!-- README.md is generated from README.Rmd. Please edit that file -->

# `{shigella}`

This repo is for the Shigella analyses

<!-- badges: start -->

[![Codecov test coverage](https://codecov.io/gh/UCD-SERG/shigella/graph/badge.svg)](https://app.codecov.io/gh/UCD-SERG/shigella) [![R-CMD-check](https://github.com/UCD-SERG/shigella/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/UCD-SERG/shigella/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

The goal of {shigella} is to implement methods for modeling longitudinal 
antibody responses to Shigella infection. The package supports hierarchical 
Bayesian modeling using JAGS, posterior summary extraction, and visualization 
of antibody decay across antigens and isotypes. It builds on and extends tools 
from the {serodynamics} and {serocalculator} packages, tailored for 
Shigella seroepidemiology.

## Installation

You can install the development version of `{shigella}` from 
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("UCD-SERG/shigella")
```