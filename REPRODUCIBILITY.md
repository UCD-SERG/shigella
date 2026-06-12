# Reproducibility notes

## Known tech debt: serodynamics internal calls

The package currently accesses the following `serodynamics` internal functions
via `:::` pending their export or replacement with public APIs:

- `ab()`
- `use_att_names()`
- `get_timeindays_var()`

`R CMD check` reports these as an accepted NOTE for this PR. The calls should
be replaced with exported `serodynamics` APIs once they are available.
