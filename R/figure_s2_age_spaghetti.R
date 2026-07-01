#' S2 Fig: raw longitudinal trajectories, children aged >= 5 years
#'
#' @param case_list_plus5 Named list of the five `>=5y` single-antigen case
#'   datasets.
#' @return A ggplot.
#' @export
figure_s2_age_spaghetti <- function(case_list_plus5) {
  figure_raw_spaghetti(case_list_plus5, numeric_y = FALSE)
}
