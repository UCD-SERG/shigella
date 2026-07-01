#' Subset case data to the combined *S. flexneri* (2a + 3a) pool
#'
#' Convenience wrapper for the cross-reactivity-informed pool used to model
#' SF2a/SF3a jointly (n = 25). Equivalent to
#' `subset_infecting_serotype(case_data, c("Sf2a", "Sf3a"))`.
#'
#' @param case_data Case data carrying a `cohort_name` column.
#' @return Subset case data with attributes preserved.
#' @export
subset_combined_flexneri <- function(case_data) {
  subset_infecting_serotype(case_data, c("Sf2a", "Sf3a"))
}
