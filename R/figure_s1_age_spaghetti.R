#' S1 Fig: raw longitudinal trajectories, children aged < 5 years
#'
#' Identical styling to main-text Figure 2 (plain log y-axis). Thin wrapper over
#' [figure_raw_spaghetti()] applied to the five `*_under5` case datasets.
#'
#' @param case_list_under5 Named list of the five `<5y` single-antigen case
#'   datasets (`IpaB`, `Sf2a`, `Sf3a`, `Sf6`, `Sonnei`), each carrying
#'   `cohort_name`.
#' @return A ggplot.
#' @export
figure_s1_age_spaghetti <- function(case_list_under5) {
  figure_raw_spaghetti(case_list_under5, numeric_y = FALSE)
}
