# Residual-based posterior predictive metrics for longitudinal antibody curves

Computes residuals between observed antibody measurements and posterior
median predictions evaluated at the observed time points. Returns
pointwise residuals or aggregated error metrics (MAE, RMSE, SSE) at
multiple summary levels.

## Usage

``` r
compute_residual_metrics(
  model,
  dataset,
  ids,
  antigen_iso,
  scale = c("original", "log"),
  summary_level = c("id_antigen", "pointwise", "antigen", "overall")
)
```

## Arguments

- model:

  A data frame of posterior draws in long format with columns:
  `Subject`, `Iso_type`, `Chain`, `Iteration`, `Parameter`, `value`.

- dataset:

  A `serodynamics` case dataset produced by
  [`serodynamics::as_case_data()`](https://rdrr.io/pkg/serodynamics/man/as_case_data.html)
  (must contain `id` and `antigen_iso` columns, and time/value
  attributes).

- ids:

  Character vector of subject IDs to include (matched against
  `dataset$id` and `model$Subject`).

- antigen_iso:

  Character scalar specifying the antigen/isotype to analyze (matched
  against `dataset$antigen_iso` and `model$Iso_type`).

- scale:

  Scale on which to compute residuals. One of `"original"` or `"log"`.
  If `"log"`, residuals are computed on the natural log scale and
  observations/predictions \\\le 0\\ are removed.

- summary_level:

  Level at which to summarize metrics. One of:

  `"pointwise"`

  :   Return pointwise residuals at each observed time point.

  `"id_antigen"`

  :   Summarize by `id` and `antigen_iso`.

  `"antigen"`

  :   Summarize by `antigen_iso` only.

  `"overall"`

  :   Single summary across all included observations.

## Value

A tibble. For `summary_level = "pointwise"`, returns per-observation
residuals. Otherwise returns MAE, RMSE, SSE, and `n_obs` at the
requested summary level.

## Details

The posterior predictive summary uses the posterior median at each
observed time point, with a 95\\ predictions. Predictions are generated
by
[`predict_posterior_at_times`](https://ucd-serg.github.io/shigella/reference/predict_posterior_at_times.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Per-ID error metrics (original scale)
m_id <- compute_residual_metrics(
  model = overall_sf2a,
  dataset = dL_clean_sf2a,
  ids = unique(dL_clean_sf2a$id),
  antigen_iso = "IgG",
  scale = "original",
  summary_level = "id_antigen"
)

# Pointwise residuals on log scale
r_pw <- compute_residual_metrics(
  model = overall_sf2a,
  dataset = dL_clean_sf2a,
  ids = "SOSAR-22008",
  antigen_iso = "IgA",
  scale = "log",
  summary_level = "pointwise"
)
} # }
```
