# Internal helpers

# Table 1 clinical variables (from the metadata workbook): sex + symptom flags.
#' @keywords internal
#' @noRd
.table1_clinical <- function(metadata) {
  metadata |>
    dplyr::transmute(
      sid = .data$sampleID,
      sex = as.character(.data$Gender),
      diarrhea_blood_mucus = factor(dplyr::case_when(
        .data$DiaBlood == 1 ~ "Yes", .data$DiaBlood == 2 ~ "No",
        TRUE ~ NA_character_
      ), levels = c("No", "Yes")),
      watery_diarrhea = factor(dplyr::case_when(
        .data$DiaWatery == 1 ~ "Yes", .data$DiaWatery == 2 ~ "No",
        TRUE ~ NA_character_
      ), levels = c("No", "Yes")),
      fever = factor(dplyr::case_when(
        .data$Fev == 1 ~ "Yes", .data$Fev == 2 ~ "No",
        TRUE ~ NA_character_
      ), levels = c("No", "Yes")),
      muac_cm = suppressWarnings(as.numeric(.data$MUAC)),
      hospital_stay_hours = suppressWarnings(as.numeric(.data$HosDur)) |>
        dplyr::na_if(88) |>
        tidyr::replace_na(0)
    ) |>
    dplyr::mutate(sex = factor(.data$sex,
      levels = c("Male", "Female", "Transgender")
    ))
}

# One row per subject: infecting serotype + age group (from the compiled sheet).
#' @keywords internal
#' @noRd
.table1_ids <- function(df_sosar) {
  df_sosar |>
    dplyr::distinct(.data$sid, .keep_all = TRUE) |>
    dplyr::mutate(
      infecting_serotype = factor(
        dplyr::case_when(
          .data$cohort_name == "Sf2a" ~ "Sf2a",
          .data$cohort_name == "sonnei" ~ "sonnei",
          .data$cohort_name == "Sf3a" ~ "Sf3a",
          .data$cohort_name == "Sf6" ~ "Sf6",
          TRUE ~ "Other"
        ),
        levels = c("Sf2a", "sonnei", "Sf3a", "Sf6", "Other")
      ),
      age_group = factor(dplyr::case_when(
        .data$age < 5 ~ "<5", .data$age >= 5 ~ "\u22655",
        TRUE ~ NA_character_
      ), levels = c("<5", "\u22655"))
    )
}

# Per-subject follow-up: visit count, follow-up days, >=4-visit flag.
#' @keywords internal
#' @noRd
.table1_followup <- function(df_sosar) {
  df_sosar |>
    dplyr::distinct(.data$sid, .data$timepoint, .keep_all = TRUE) |>
    dplyr::group_by(.data$sid) |>
    dplyr::summarise(
      n_visits = dplyr::n(),
      followup_days = max(.data[["Actual day"]], na.rm = TRUE) + 2, # +2: days since symptom onset to day-0 blood draw # nolint: line_length_linter.
      .groups = "drop"
    ) |>
    dplyr::mutate(ge4_visits = factor(dplyr::if_else(.data$n_visits >= 4, "\u22654", "<4"), # nolint: line_length_linter.
      levels = c("<4", "\u22654")
    ))
}

# Main function

#' Assemble the participant-level data frame underlying Table 1
#'
#' Joins per-subject ids (`.table1_ids()`), follow-up (`.table1_followup()`),
#' clinical variables (`.table1_clinical()`) and pre-presentation diarrhea
#' duration onto the SOSAR cohort, then derives `completed_followup`.
#'
#' @param compiled `Compiled` sheet (MFI + `cohort_name` + `age`).
#' @param metadata `Metadata` sheet (sex, clinical signs, MUAC, hospital stay).
#' @param durdia Diarrhea-duration sheet (`CaseID`, `DurDia_hours`).
#' @return A participant-level tibble ready for [table1_study_population()].
#' @export
build_table1_data <- function(compiled, metadata, durdia) {
  df_sosar <- dplyr::filter(compiled, .data$study_name == "SOSAR")

  df_durdia <- durdia |>
    dplyr::transmute(
      sid = .data$CaseID,
      durdia_hours = suppressWarnings(as.numeric(.data$DurDia_hours))
    ) # nolint: line_length_linter.

  .table1_ids(df_sosar) |>
    dplyr::left_join(.table1_followup(df_sosar), by = "sid") |>
    dplyr::left_join(.table1_clinical(metadata), by = "sid") |>
    dplyr::left_join(df_durdia, by = "sid") |>
    dplyr::mutate(
      sex = droplevels(.data$sex),
      completed_followup = factor(dplyr::case_when(
        .data$n_visits >= 4 & .data$followup_days >= 90 ~ "Yes",
        TRUE ~ "No"
      ), levels = c("Yes", "No"))
    )
}
