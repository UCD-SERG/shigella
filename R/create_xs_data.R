#' Create Cross-Sectional Data for Specific Regions and Antigens
#'
#' Filters and transforms Shigella study data to create a cross-sectional dataset
#' for specific geographic regions and antigen-isotype combinations, including
#' age categorization.
#'
#' @param df A data frame containing Shigella study data with columns:
#'   `sid` (subject ID), `site_name`, `study_name`, `age`, `isotype_name`,
#'   and antigen-specific measurement columns.
#' @param filter_countries Character vector of country/site names to include
#'   (e.g., `c("MA USA", "Ghana", "Niger")`).
#' @param filter_antigen_iso Character vector of antigen-isotype combinations
#'   to include (e.g., `c("IgG", "IgA")`).
#' @param value_col Unquoted column name containing the antigen measurement
#'   values to extract (e.g., `n_ipab_MFI`, `n_sf2a_MFI`).
#'
#' @return A data frame with cross-sectional data containing columns:
#'   `id`, `Country`, `study`, `age`, `antigen_iso`, `value`, and `ageCat`
#'   (age category: "<5", "5-15", or "16+").
#'
#' @examples
#' \dontrun{
#' # Create mock Shigella data
#' mock_df <- data.frame(
#'   sid = 1:100,
#'   site_name = sample(c("MA USA", "Ghana", "Niger"), 100, replace = TRUE),
#'   study_name = rep("SOSAR", 100),
#'   age = runif(100, 0.5, 50),
#'   isotype_name = sample(c("IgG", "IgA"), 100, replace = TRUE),
#'   n_ipab_MFI = rlnorm(100, log(1000), 1)
#' )
#'
#' # Extract cross-sectional data for USA IgG
#' xs_data <- create_xs_data(
#'   df = mock_df,
#'   filter_countries = c("MA USA"),
#'   filter_antigen_iso = c("IgG"),
#'   value_col = n_ipab_MFI
#' )
#' }
#'
#' @importFrom dplyr mutate filter select %>%
#' @importFrom rlang {{ }}
#' @export
# create function for generating specific region data
create_xs_data <- function(df, filter_countries, filter_antigen_iso, value_col) {
  df %>%
    # First, create/rename columns
    mutate(
      id = sid,
      Country = site_name,
      study = study_name,
      age = age,
      antigen_iso = factor(isotype_name),
      value = {{ value_col }
      } # value_col is specified by the user, e.g. n_ipab_MFI
    ) %>%
    # Then filter by the desired catchment(s) and antigen iso(s)
    filter(
      Country %in% filter_countries,
      antigen_iso %in% filter_antigen_iso
    ) %>%
    # Create age categories
    mutate(
      ageCat = factor(case_when(
        age < 5 ~ "<5",
        age >= 5 & age <= 15 ~ "5-15",
        age > 15 ~ "16+"
      ))
    ) %>%
    # Optionally select only the needed columns
    select(id, Country, study, age, antigen_iso, value, ageCat)
}