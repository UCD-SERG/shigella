# -----------------------------------------------------------------------------
# Internal: apply y0/y1/t1/alpha/rho sub/superscript headers to a flextable
# (only to the columns that are present). Reused by Tables 2 and 3.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.add_kinetic_headers <- function(ft) {
  cols <- ft$col_keys
  comp <- function(ft, j, value) {
    if (j %in% cols) flextable::compose(ft, part = "header", j = j, value = value) else ft # nolint: line_length_linter.
  }
  ft |>
    comp("y0",    flextable::as_paragraph(flextable::as_chunk("y"), flextable::as_sub("0"))) |> # nolint: line_length_linter.
    comp("y1",    flextable::as_paragraph(flextable::as_chunk("y"), flextable::as_sub("1"))) |> # nolint: line_length_linter.
    comp("t1",    flextable::as_paragraph(flextable::as_chunk("t"), flextable::as_sub("1"))) |> # nolint: line_length_linter.
    comp("alpha", flextable::as_paragraph(flextable::as_chunk("\u03B1"),
                                          flextable::as_chunk(" (year"),
                                          flextable::as_sup("-1"),
                                          flextable::as_chunk(")"))) |>
    comp("rho",   flextable::as_paragraph(flextable::as_chunk("\u03C1")))
}
