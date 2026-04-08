# Compare serotype-specific vs overall models using residual metrics

Computes per-ID residual metrics for two models on the intersection of
IDs present in both datasets, then reports absolute and percent
differences.

## Usage

``` r
make_model_comparison_table(
  model_serospecific,
  data_serospecific,
  model_overall,
  data_overall,
  antigen_iso,
  scale = c("original", "log"),
  tie_tol = 1e-08
)
```

## Arguments

- model_serospecific:

  Posterior draws (long format) for the serotype-specific model.

- data_serospecific:

  Case dataset used for the serotype-specific model.

- model_overall:

  Posterior draws (long format) for the overall model.

- data_overall:

  Case dataset used for the overall model.

- antigen_iso:

  Character scalar antigen/isotype label.

- scale:

  "original" or "log".

- tie_tol:

  Numeric tolerance to declare a tie.

## Value

A tibble with per-ID MAE/RMSE for each model and deltas/winner labels.
