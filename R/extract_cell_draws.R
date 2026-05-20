# Helper: extract all chain draws for one parameter cell `pname[subj, iso]`.
# Parses the column name for indices, slices to the requested chains, and
# returns a list of one tibble per chain (ready for dplyr::bind_rows).
#' @keywords internal
#' @noRd
.extract_draws_for_cell <- function(col_name, draws_df, pname,
                                    ids, antigens, stratification, n_chain) {
  m        <- regmatches(col_name, regexec("\\[(\\d+),(\\d+)\\]", 
                                           col_name))[[1]]
  subj_idx <- as.integer(m[2])
  iso_idx  <- as.integer(m[3])
  sub_df   <- draws_df[, c(".chain", ".iteration", col_name)]

  lapply(seq_len(n_chain), function(ch) {
    chain_data <- sub_df[sub_df$.chain == ch, ]
    tibble::tibble(
      Iteration      = chain_data$.iteration,
      Chain          = ch,
      Parameter      = pname,
      Iso_type       = antigens[iso_idx],
      Stratification = stratification,
      Subject        = ids[subj_idx],
      value          = chain_data[[col_name]]
    )
  })
}
