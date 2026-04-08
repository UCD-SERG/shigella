# Extract "newperson" parameter draws and summarize for Table 2

Extract "newperson" parameter draws and summarize for Table 2

## Usage

``` r
prep_newperson_params(draws_long, antigen_label)
```

## Arguments

- draws_long:

  Posterior draws in long format with columns: Subject, Iso_type, Chain,
  Iteration, Parameter, value.

- antigen_label:

  Character scalar antigen name to attach.

## Value

A tibble of newperson draws in wide format with columns: Iteration,
Chain, antigen, Iso_type, y0, y1, t1, alpha, rho.
