#' Process Shigella Longitudinal Data
#'
#' Filters and transforms raw Shigella study data for a specific study and antigen,
#' preparing it for downstream analysis and visualization.
#'
#' @param data A data frame containing raw Shigella study data with columns:
#'   `study_name`, `isotype_name`, `sid` (subject ID), `timepoint`,
#'   `Actual day`, and antigen-specific columns.
#' @param study_filter Character string specifying the study name to filter by
#'   (e.g., "SOSAR", "GEMS").
#' @param antigen Unquoted column name of the antigen measurement to extract
#'   (e.g., `ipab_MFI`, `sf2a_MFI`).
#'
#' @return A data frame with processed longitudinal data containing columns:
#'   `index_id`, `antigen_iso`, `visit`, `timeindays`, `result`, and `visit_num`.
#'   Rows with missing `timeindays` are removed.
#'
#' @examples
#' \dontrun{
#' # Create mock data
#' mock_data <- data.frame(
#'   study_name = rep("SOSAR", 10),
#'   isotype_name = rep(c("IgG", "IgA"), each = 5),
#'   sid = rep(1:5, 2),
#'   timepoint = rep(c("V1", "V2"), 5),
#'   `Actual day` = c(0, 14, 28, 42, 56, 0, 14, 28, 42, 56),
#'   ipab_MFI = rnorm(10, 1000, 200)
#' )
#'
#' # Process the data
#' processed <- process_shigella_data(
#'   data = mock_data,
#'   study_filter = "SOSAR",
#'   antigen = ipab_MFI
#' )
#' }
#'
#' @importFrom dplyr filter select mutate group_by arrange ungroup %>%
#' @importFrom rlang ensym !!
#' @export
process_shigella_data <- function(data, study_filter, antigen) {
  # Filter the data for the specific study
  filtered_data <- data %>%
    filter(study_name == study_filter)
  
  # Capture the column name of the antigen
  antigen_col <- ensym(antigen)
  
  # Manipulate and restructure the data
  processed_data <- filtered_data %>%
    select(isotype_name, sid, timepoint, `Actual day`, !!antigen_col) %>%
    mutate(
      index_id = sid,
      antigen_iso = isotype_name,
      visit = timepoint,
      timeindays = `Actual day`,
      result = !!antigen_col
    ) %>%
    group_by(index_id, antigen_iso) %>%
    arrange(visit) %>%
    mutate(visit_num = rank(visit, ties.method = "first")) %>%
    ungroup() %>%
    # Remove rows with NA in timeindays
    filter(!is.na(timeindays))
  
  return(processed_data)
}