#' Table 3: age-stratified IpaB population-level kinetic parameters
#'
#' @param model_under5,model_plus5 Fitted IpaB `sr_model`s for each age group.
#' @return A flextable (biomarker outer, age group nested).
#' @export
table3_age_stratified <- function(model_under5, model_plus5) {
  age_summary <- dplyr::bind_rows(
    dplyr::mutate(extract_mu_draws(model_under5, "IpaB"), age_group = "<5"),
    dplyr::mutate(extract_mu_draws(model_plus5,  "IpaB"), age_group = "\u22655")
  ) |>
    dplyr::mutate(Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA"))) |>
    add_natural_scale(alpha_unit = "per_year") |>
    summarise_pop_params(group_vars = c("age_group", "antigen", "Iso_type"))

  age_print <- age_summary |>
    dplyr::mutate(Biomarker = paste0("IpaB\u2013", .data$Iso_type)) |>
    dplyr::arrange(.data$Iso_type, .data$age_group) |>
    dplyr::mutate(
      `Age group` = .data$age_group,
      y0    = fmt_mci(.data$y0_med, .data$y0_lo, .data$y0_hi, digits = 2),
      y1    = fmt_mci(.data$y1_med, .data$y1_lo, .data$y1_hi, digits = 2),
      t1    = fmt_mci(.data$t1_med, .data$t1_lo, .data$t1_hi, digits = 1),
      alpha = fmt_mci(.data$alpha_med, .data$alpha_lo, .data$alpha_hi, digits = 8),
      rho   = fmt_mci(.data$rho_med, .data$rho_lo, .data$rho_hi, digits = 2)) |>
    dplyr::select("Biomarker", "Age group", "y0", "y1", "t1", "alpha", "rho")

  ft <- age_print |>
    flextable::flextable() |>
    flextable::theme_booktabs() |>
    flextable::fontsize(size = 7, part = "body") |>
    flextable::fontsize(size = 8, part = "header") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = 1:2, align = "left", part = "body") |>
    flextable::merge_v(j = "Biomarker") |>
    flextable::set_table_properties(layout = "autofit", width = 1) |>
    flextable::set_caption(paste0(
      "Age-stratified population-level antibody kinetic parameters for IpaB ",
      "from the overall hierarchical model.")) |>
    flextable::add_footer_lines(values = c(
      "Values are posterior medians (95% credible intervals) for population mean (mu.par) draws.",
      "Age stratification shown only for IpaB (conserved antigen); O-antigen responses require serotype conditioning, producing subgroups too small for stable estimation.")) |>
    flextable::fontsize(size = 6, part = "footer") |>
    flextable::italic(part = "footer") |>
    flextable::color(color = "gray40", part = "footer")

  .add_kinetic_headers(ft) # nolint: object_usage_linter
}
