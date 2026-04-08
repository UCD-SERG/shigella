# Mock case data for testing

A mock longitudinal antibody dataset compatible with
[`serodynamics::as_case_data()`](https://rdrr.io/pkg/serodynamics/man/as_case_data.html),
for testing functions like
[`compute_residual_metrics`](https://ucd-serg.github.io/shigella/reference/compute_residual_metrics.md).

## Usage

``` r
mock_case_data
```

## Format

A data frame with class `c("case_data", "data.frame")` and columns:

- id:

  Character. Subject ID

- antigen_iso:

  Character. Isotype (e.g., "IgG", "IgA")

- timepoint:

  Numeric. Time in days since infection

- value:

  Numeric. Antibody measurement

## Details

This is synthetic data generated for testing and examples. Real Shigella
case data will be added separately.

The dataset has attributes:

- `attr(mock_case_data, "timeindays") = "timepoint"`

- `attr(mock_case_data, "value_var") = "value"`

## Examples

``` r
head(mock_case_data)
#>            id antigen_iso timepoint    value
#> 1 SOSAR-22008         IgG         0 241.3326
#> 2 SOSAR-22015         IgG         0 152.5996
#> 3 SOSAR-22020         IgG         0 272.4704
#> 4 SOSAR-22025         IgG         0 327.1269
#> 5 SOSAR-22008         IgA         0 255.5663
#> 6 SOSAR-22015         IgA         0 277.2024
table(mock_case_data$id)
#> 
#> SOSAR-22008 SOSAR-22015 SOSAR-22020 SOSAR-22025 
#>          12          12          12          12 
```
