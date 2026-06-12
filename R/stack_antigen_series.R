#' Stack single-antigen case data into one long frame for spaghetti plots
#'
#' @param case_list Named list of single-antigen case-data objects (carrying
#'   `cohort_name`), names giving the antigen labels
#'   (`IpaB`, `Sf2a`, `Sf3a`, `Sf6`, `Sonnei`).
#' @return A tibble with `antigen`, `isotype_name`, `infecting_serotype`
#'   factors.
#' @export
stack_antigen_series <- function(case_list) {
  purrr::imap(case_list, ~ dplyr::mutate(.x, antigen = .y)) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      antigen      = factor(.data$antigen,
                            levels = c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei")), # nolint: line_length_linter.
      isotype_name = factor(.data$isotype_name, levels = c("IgG", "IgA")),
      infecting_serotype = dplyr::case_when(
        .data$cohort_name == "Sf2a"   ~ "Sf2a",
        .data$cohort_name == "sonnei" ~ "sonnei",
        .data$cohort_name == "Sf3a"   ~ "Sf3a",
        .data$cohort_name == "Sf6"    ~ "Sf6",
        TRUE                          ~ "Other"
      ),
      infecting_serotype = factor(.data$infecting_serotype,
                                  levels = c("Sf2a", "sonnei", "Sf3a", "Sf6", "Other")) # nolint: line_length_linter.
    )
}
