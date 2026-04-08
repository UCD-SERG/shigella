# Mock posterior draws for testing

A mock dataset of posterior parameter draws in long format, mimicking
the structure expected by functions like
[`compute_residual_metrics`](https://ucd-serg.github.io/shigella/reference/compute_residual_metrics.md)
and
[`predict_posterior_at_times`](https://ucd-serg.github.io/shigella/reference/predict_posterior_at_times.md).

## Usage

``` r
mock_posterior_draws
```

## Format

A data frame with columns:

- Subject:

  Character. Subject ID (e.g., "newperson", "SOSAR-22008")

- Iso_type:

  Character. Isotype (e.g., "IgG", "IgA")

- Chain:

  Integer. MCMC chain number

- Iteration:

  Integer. MCMC iteration number

- Parameter:

  Character. Parameter name (y0, y1, t1, alpha, shape)

- value:

  Numeric. Parameter value

## Details

This is synthetic data generated for testing and examples. Real Shigella
posterior draws will be added to the package separately.

Parameters represent:

- `y0`: Baseline antibody level

- `y1`: Peak antibody level

- `t1`: Time to peak (days)

- `alpha`: Decay rate parameter

- `shape`: Decay shape parameter (rho)

## Examples

``` r
head(mock_posterior_draws)
#>       Subject Iso_type Chain Iteration Parameter    value
#> 1   newperson      IgG     1         1        y0 434.7770
#> 2 SOSAR-22008      IgG     1         1        y0 228.3470
#> 3 SOSAR-22015      IgG     1         1        y0 372.1453
#> 4   newperson      IgA     1         1        y0 379.2692
#> 5 SOSAR-22008      IgA     1         1        y0 282.8037
#> 6 SOSAR-22015      IgA     1         1        y0 380.5681
table(mock_posterior_draws$Parameter)
#> 
#> alpha shape    t1    y0    y1 
#>  1800  1800  1800  1800  1800 
```
