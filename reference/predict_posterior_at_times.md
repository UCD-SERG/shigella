# Posterior predictions at specified times for given subjects and antigen/isotype

Generates draw-level posterior predictions of the antibody trajectory at
user-specified time points, for one or more subjects and a selected
antigen/isotype. This is a low-level helper used by residual-based
posterior predictive diagnostics.

## Usage

``` r
predict_posterior_at_times(model, ids, antigen_iso, times)
```

## Arguments

- model:

  A data frame of posterior draws in long format with columns:
  `Subject`, `Iso_type`, `Chain`, `Iteration`, `Parameter`, `value`.

- ids:

  Character vector of subject IDs to include (matched against
  `Subject`).

- antigen_iso:

  Character scalar specifying the antigen/isotype to include (matched
  against `Iso_type`).

- times:

  Numeric vector of time points (days) at which to evaluate predictions.

## Value

A tibble with one row per (posterior draw \\\times\\ time \\\times\\
subject), including the evaluated prediction `res`. Output includes at
least:

- id:

  Subject ID (character).

- t:

  Time (days) at which prediction was evaluated.

- Chain:

  MCMC chain index (if present in `model`).

- Iteration:

  MCMC iteration index (if present in `model`).

- sample_id:

  Row index for the draw (added if missing).

- y0, y1, t1, alpha, shape:

  Model parameters (wide).

- res:

  Predicted antibody level at time `t`.

## Details

This function pivots posterior draws to wide format (parameters as
columns), expands them over `times`, and evaluates the antibody curve
via an internal implementation of the antibody kinetics model using
parameters `y0`, `y1`, `t1`, `alpha`, and `shape`.

## See also

[`compute_residual_metrics`](https://ucd-serg.github.io/shigella/reference/compute_residual_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
preds <- predict_posterior_at_times(
  model = overall_sf2a,
  ids = "newperson",
  antigen_iso = "IgG",
  times = c(0, 30, 90, 180)
)
} # }
```
