#' Table 1: cohort characteristics (gtsummary -> flextable)
#'
#' @param tab1_dat Output of [build_table1_data()].
#' @return A flextable.
#' @export
table1_study_population <- function(tab1_dat) {
  tab1_gts <- tab1_dat |>
    dplyr::select(
      "infecting_serotype", "age", "age_group", "sex",
      "diarrhea_blood_mucus", "watery_diarrhea", "fever", "muac_cm",
      "hospital_stay_hours", "durdia_hours", "n_visits", "followup_days",
      "completed_followup") |>
    gtsummary::tbl_summary(
      by = "infecting_serotype",
      type = list(sex ~ "dichotomous"),
      value = list(sex ~ "Male"),
      statistic = list(
        gtsummary::all_continuous() ~ "{median} ({p25}\u2013{p75})",
        hospital_stay_hours ~ "{median} ({min}\u2013{max})",
        followup_days ~ "{median} ({min}\u2013{max})",
        gtsummary::all_categorical() ~ "{n} ({p}%)"),
      digits = gtsummary::all_continuous() ~ 1,
      missing = "no",
      label = list(
        age ~ "Age at enrollment, years",
        age_group ~ "Age group, years",
        sex ~ "Male sex, n (%)",
        diarrhea_blood_mucus ~ "Diarrhea with blood or mucus",
        watery_diarrhea ~ "Watery diarrhea",
        fever ~ "Fever",
        muac_cm ~ "MUAC (cm)",
        hospital_stay_hours ~ "Duration of hospital stay (hours)",
        durdia_hours ~ "Diarrhea duration prior to presentation (hours)",
        n_visits ~ "Samples per participant",
        followup_days ~ "Follow-up duration, days",
        completed_followup ~ "Completed follow-up (\u22654 visits and \u226590 days)")) |>
    gtsummary::add_overall(last = FALSE, col_label = "**Overall**") |>
    gtsummary::modify_header(label ~ "**Characteristic**") |>
    gtsummary::bold_labels() |>
    gtsummary::modify_spanning_header(gtsummary::all_stat_cols() ~ NA_character_) |>
    gtsummary::modify_header(stat_0 ~ paste0("**Overall (N = ", nrow(tab1_dat), ")**"))

  tab1_gts |>
    gtsummary::as_flex_table() |>
    flextable::bold(part = "header") |>
    flextable::fontsize(size = 7, part = "all") |>
    flextable::fontsize(size = 8, part = "header") |>
    flextable::set_table_properties(layout = "autofit", width = 1) |>
    flextable::set_caption("Characteristics of the longitudinal Shigella cohort.") |>
    flextable::add_footer_lines(values = c(
      "Continuous variables: median (Q1\u2013Q3) except hospital stay and follow-up: median (min\u2013max).",
      "Follow-up duration (days) computed as max(Actual day) + 2 to account for enrollment on day 2 after presentation.",
      "\"Other\" serotypes include S. boydii (n=2), S. dysenteriae (n=1), S. flexneri 1c (n=1), untypeable S. flexneri (n=2), and RLDT-positive Shigella without culture confirmation (n=1).",
      "Not all 48 participants completed full follow-up; the \"Other\" group had shorter median follow-up (32 days) due to loss to follow-up."))
}
