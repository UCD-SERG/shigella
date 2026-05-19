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
  out_list <- list()
  row_counter <- 1L

  for (pname in param_names) {
    # Find columns matching pname[i,k]
    matching_cols <- grep(paste0("^", pname, "\\["), colnames(draws_df),
                          value = TRUE)
    if (length(matching_cols) != N * K) {
      cli::cli_abort(
        "Expected {N * K} {.var {pname}} draws; got {length(matching_cols)}."
      )
    }

    for (col_name in matching_cols) {
      m <- regmatches(col_name, regexec("\\[(\\d+),(\\d+)\\]", col_name))[[1]]
      subj_idx <- as.integer(m[2])
      iso_idx  <- as.integer(m[3])

      # Extract draws for this parameter index
      sub_df <- draws_df[, c(".chain", ".iteration", col_name)]

      for (ch in seq_len(n_chain)) {
        chain_data <- sub_df[sub_df$.chain == ch, ]
        out_list[[row_counter]] <- tibble::tibble(
          Iteration      = chain_data$.iteration,
          Chain          = ch,
          Parameter      = pname,
          Iso_type       = antigens[iso_idx],
          Stratification = stratification,
          Subject        = ids[subj_idx],
          value          = chain_data[[col_name]]
        )
        row_counter <- row_counter + 1L
      }
    }
  }

  dplyr::bind_rows(out_list)
}
