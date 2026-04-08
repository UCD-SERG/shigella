# Getting Started with shigella

``` r
library(shigella)
```

## Overview

The `shigella` package provides tools for analyzing longitudinal
antibody kinetics data from Shigella infection studies. This vignette
introduces the main functionality using mock data.

**Note:** This package currently uses mock data for testing and
examples. Real Shigella datasets will be added in a future release.

## Key Functions

### Data Processing

The
[`process_shigella_data()`](https://ucd-serg.github.io/shigella/reference/process_shigella_data.md)
function reshapes raw longitudinal data into a format compatible with
the `serodynamics` package:

``` r
# Example with real data (not run)
dat_long <- process_shigella_data(
  data = raw_data,
  study_filter = "SOSAR",
  antigen = n_ipab_MFI
)

# Convert to case_data format
dL <- serodynamics::as_case_data(
  dat_long,
  id_var = "index_id",
  biomarker_var = "antigen_iso",
  time_in_days = "timeindays",
  value_var = "result"
)
```

### Model Evaluation

The package provides functions to evaluate and compare longitudinal
antibody models:

#### Computing Residual Metrics

``` r
# Compute residuals at individual level
metrics_id <- compute_residual_metrics(
  model = posterior_draws,
  dataset = case_data,
  ids = unique(case_data$id),
  antigen_iso = "IgG",
  scale = "original",
  summary_level = "id_antigen"
)

# Overall summary metrics
metrics_overall <- compute_residual_metrics(
  model = posterior_draws,
  dataset = case_data,
  ids = unique(case_data$id),
  antigen_iso = "IgG",
  scale = "original",
  summary_level = "overall"
)
```

#### Model Comparison

Compare model fit metrics between different modeling approaches:

``` r
comparison <- model_comparison_table(
  metrics_overall = metrics_model1,
  metrics_pointwise = metrics_model2,
  model_overall_label = "Population Model",
  model_pointwise_label = "Individual Model"
)
```

### Visualization

The
[`fig2_overall_newperson()`](https://ucd-serg.github.io/shigella/reference/fig2_overall_newperson.md)
function creates trajectory plots for population-level (“newperson”)
predictions:

``` r
fig <- fig2_overall_newperson(
  overall_models = list(
    IpaB = model_ipab,
    Sf2a = model_sf2a
  ),
  osps = c("IpaB", "Sf2a"),
  isotypes = c("IgG", "IgA"),
  t_grid = seq(0, 210, by = 5),
  log_y = TRUE
)
```

## Mock Data

The package includes two mock datasets for testing:

- `mock_posterior_draws`: Posterior parameter draws in long format
- `mock_case_data`: Longitudinal antibody measurements

These datasets have the same structure as real data but contain
synthetic values.

``` r
# Load mock data
data("mock_posterior_draws")
data("mock_case_data")

# View structure
head(mock_posterior_draws)
head(mock_case_data)
```

## Next Steps

For more detailed examples and analysis workflows, see:

- Package documentation:
  [`?shigella`](https://ucd-serg.github.io/shigella/reference/shigella-package.md)
- Function help:
  [`?compute_residual_metrics`](https://ucd-serg.github.io/shigella/reference/compute_residual_metrics.md),
  [`?process_shigella_data`](https://ucd-serg.github.io/shigella/reference/process_shigella_data.md)
- GitHub repository: <https://github.com/UCD-SERG/shigella>
