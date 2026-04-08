# Format median and credible interval as a single string

Format median and credible interval as a single string

## Usage

``` r
fmt_mci(med, lo, hi, digits = 2, sci = FALSE)
```

## Arguments

- med:

  Median value.

- lo:

  Lower bound.

- hi:

  Upper bound.

- digits:

  Number of digits.

- sci:

  Logical; if TRUE uses scientific notation.

## Value

A string like "1.23 (0.50–2.00)".
