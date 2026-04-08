# Compare summary residual metrics between two model outputs

Builds a compact comparison table from pre-computed summary metrics for
an overall model and a pointwise model.

## Usage

``` r
model_comparison_table(
  metrics_overall,
  metrics_pointwise,
  model_overall_label = "Overall Model",
  model_pointwise_label = "Pointwise Model"
)
```

## Arguments

- metrics_overall:

  Data frame with one row containing at least `MAE`, `RMSE`, `SSE`, and
  `n_obs` for the overall model.

- metrics_pointwise:

  Data frame with one row containing at least `MAE`, `RMSE`, `SSE`, and
  `n_obs` for the pointwise model.

- model_overall_label:

  Label used for the overall model row.

- model_pointwise_label:

  Label used for the pointwise model row.

## Value

A tibble with one row per model and comparison columns: `delta_MAE`,
`delta_RMSE`, `pct_improve_MAE`, `pct_improve_RMSE`.
