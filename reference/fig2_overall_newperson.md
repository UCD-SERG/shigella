# Summarize population-level ("newperson") antibody trajectories from overall models

Computes median and credible interval bands of the antibody trajectory
for a hypothetical new individual drawn from the population distribution
("newperson").

## Usage

``` r
fig2_overall_newperson(
  overall_models,
  osps = c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei"),
  ids = "newperson",
  isotypes = c("IgG", "IgA"),
  t_grid = seq(0, 210, by = 5),
  cred = 0.95,
  log_y = TRUE,
  xlim = c(0, 210),
  ylab = "Normalized MFI",
  line_color = "#1f77b4",
  ribbon_alpha = 0.2,
  facet_scales = "fixed",
  return_data = FALSE
)
```

## Arguments

- overall_models:

  Named list of posterior draws in long format (one per antigen), with
  columns at least: Subject, Iso_type, Chain, Iteration, Parameter,
  value.

- osps:

  Character vector of antigen names (must match names in
  `overall_models`).

- ids:

  Character vector of Subject IDs to include (default: "newperson").

- isotypes:

  Character vector of isotypes (default: c("IgG", "IgA")).

- t_grid:

  Numeric vector of time points (days) to evaluate.

- cred:

  Credible level (default 0.95).

- log_y:

  Logical; if TRUE, applies log10 scale to y when returning plot.

- xlim:

  Optional numeric length-2 vector for x-axis limits when returning
  plot.

- ylab:

  Y-axis label for plot.

- line_color:

  Line color for plot.

- ribbon_alpha:

  Alpha for credible ribbon in plot.

- facet_scales:

  Passed to ggplot facet scales.

- return_data:

  If TRUE, returns a list with `plot` and `data`.

## Value

By default, a ggplot. If `return_data = TRUE`, returns list(plot = p,
data = df).

## Details

Requires that each model draw can be pivoted to wide parameters
including: y0, y1, t1, alpha, shape.
