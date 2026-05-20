# Helper: extract draws for all param_names into a single tibble.
#' @keywords internal
#' @noRd
.extract_param_draws <- function(param_names,
                                 draws_df,
                                 N,
                                 K,
                                 ids,
                                 antigens,
                                 stratification,
                                 n_chain) {
  draws_per_param <- lapply(param_names, function(pname) {
    matching_cols <- grep(paste0("^", pname, "\\["), colnames(draws_df),
                          value = TRUE)
    if (length(matching_cols) != N * K) {
      cli::cli_abort(
        "Expected {N * K} {.var {pname}} draws; got {length(matching_cols)}."
      )
    }
    draws_per_col <- lapply(
      matching_cols, .extract_draws_for_cell,
      draws_df       = draws_df,
      pname          = pname,
      ids            = ids,
      antigens       = antigens,
      stratification = stratification,
      n_chain        = n_chain
    )
    dplyr::bind_rows(unlist(draws_per_col, recursive = FALSE))
  })

  dplyr::bind_rows(draws_per_param)
}
