# Reshape Shigella longitudinal data for serodynamics workflows

Filters a dataset to a given study and antigen column, standardizes
column names, and returns a visit-ordered long dataset suitable for
conversion to
[`serodynamics::as_case_data()`](https://rdrr.io/pkg/serodynamics/man/as_case_data.html).

## Usage

``` r
process_shigella_data(data, study_filter, antigen)
```

## Arguments

- data:

  A data frame containing longitudinal measurements.

- study_filter:

  Character scalar. Value of `study_name` to keep (e.g. "SOSAR").

- antigen:

  Unquoted column name for the antigen measurement (e.g. n_ipab_MFI).

## Value

A tibble with standardized columns:

- index_id:

  Participant ID (copied from `sid`).

- antigen_iso:

  Isotype label (copied from `isotype_name`).

- visit:

  Visit label (copied from `timepoint`).

- timeindays:

  Time since infection (copied from `Actual day`).

- result:

  Antibody measurement (from `antigen`).

## Examples

``` r
if (FALSE) { # \dontrun{
dat_long <- process_shigella_data(df, "SOSAR", n_ipab_MFI)
dL <- serodynamics::as_case_data(dat_long,
  id_var = "index_id", biomarker_var = "antigen_iso",
  time_in_days = "timeindays", value_var = "result"
)
} # }
```
