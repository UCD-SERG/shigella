# ============================================================
# data-raw/mock_data.R
#
# Generates lightweight mock versions of all data objects used
# by the shigella package. These mock objects have identical
# STRUCTURE to the real JAGS model outputs and case-data objects
# but contain only 5 subjects, 5 timepoints, and 20 MCMC
# iterations (instead of 48 subjects × 60,000 iterations).
#
# Run this script from the package root:
#   source("data-raw/mock_data.R")
#
# It will populate data/ with .rda files that allow:
#   - R CMD check to pass
#   - Tests to run
#   - Examples and vignettes to execute
# ============================================================

library(tibble)
library(dplyr)

# ── Configuration ──────────────────────────────────────────────────────────

set.seed(42)

# Mock subjects: 3 Sf2a, 1 Sonnei, 1 Sf3a
mock_subjects <- tibble::tibble(

  sid         = paste0("MOCK-", sprintf("%03d", 1:5)),
  cohort_name = c("Sf2a", "Sf2a", "sonnei", "Sf3a", "Sf2a"),
  age         = c(3, 8, 2, 15, 25)
)

# Timepoints (days since symptom onset)
timepoints_nominal <- c(2, 7, 30, 90, 180)

# MCMC dimensions (tiny for mock)
n_iter   <- 20L
n_chains <- 2L

# Isotypes
isotypes <- c("IgG", "IgA")

# Parameters (on log scale, as stored by JAGS)
# Realistic values based on IpaB from the manuscript
param_names <- c("y0", "y1", "t1", "alpha", "shape")

# Population-level means (log scale) for different antigens/isotypes
# These produce realistic natural-scale values when exponentiated
pop_means <- list(
  IpaB_IgG   = c(y0 = 5.0, y1 = 10.0, t1 = -0.5, alpha = -15.0, shape = 0.3),
  IpaB_IgA   = c(y0 = 3.5, y1 =  9.0, t1 = -0.8, alpha = -12.0, shape = 0.2),
  Sf2a_IgG   = c(y0 = 6.0, y1 =  9.5, t1 = -0.3, alpha = -14.0, shape = 0.4),
  Sf2a_IgA   = c(y0 = 4.0, y1 =  8.5, t1 = -0.6, alpha = -11.0, shape = 0.3),
  Sf3a_IgG   = c(y0 = 5.5, y1 =  9.0, t1 = -0.4, alpha = -13.5, shape = 0.35),
  Sf3a_IgA   = c(y0 = 3.8, y1 =  8.0, t1 = -0.7, alpha = -11.5, shape = 0.25),
  Sf6_IgG    = c(y0 = 4.5, y1 =  7.5, t1 = -0.5, alpha = -14.5, shape = 0.3),
  Sf6_IgA    = c(y0 = 3.0, y1 =  6.5, t1 = -0.8, alpha = -12.5, shape = 0.2),
  Sonnei_IgG = c(y0 = 5.5, y1 = 10.5, t1 = -0.3, alpha = -13.0, shape = 0.5),
  Sonnei_IgA = c(y0 = 4.0, y1 =  9.5, t1 = -0.5, alpha = -10.0, shape = 0.4)
)


# ── Helper: create sr_model object ────────────────────────────────────────

#' @param subjects Character vector of subject IDs to include
#' @param antigen_means Named list of pop means (only IgG/IgA entries used)
#' @param antigen_key Prefix for selecting from pop_means (e.g., "IpaB")
make_mock_model <- function(subjects, antigen_key, all_pop_means = pop_means) {

  igg_means <- all_pop_means[[paste0(antigen_key, "_IgG")]]
  iga_means <- all_pop_means[[paste0(antigen_key, "_IgA")]]

  # --- Individual-level draws ---
  indiv_rows <- list()
  for (sid in subjects) {
    for (iso in isotypes) {
      mu <- if (iso == "IgG") igg_means else iga_means
      for (ch in seq_len(n_chains)) {
        for (it in seq_len(n_iter)) {
          # Draw individual params with small noise around population mean
          vals <- mu + rnorm(5, 0, 0.3)
          names(vals) <- param_names
          for (p in param_names) {
            indiv_rows[[length(indiv_rows) + 1]] <- tibble::tibble(
              Iteration      = it,
              Chain          = ch,
              Parameter      = p,
              Iso_type       = iso,
              Stratification = "None",
              Subject        = sid,
              value          = vals[p]
            )
          }
        }
      }
    }
  }
  indiv_df <- dplyr::bind_rows(indiv_rows)

  # --- Population-level draws (mu.par) ---
  pop_rows <- list()
  for (iso in isotypes) {
    mu <- if (iso == "IgG") igg_means else iga_means
    for (ch in seq_len(n_chains)) {
      for (it in seq_len(n_iter)) {
        vals <- mu + rnorm(5, 0, 0.15)
        names(vals) <- param_names
        for (p in param_names) {
          pop_rows[[length(pop_rows) + 1]] <- tibble::tibble(
            Iteration            = it,
            Chain                = ch,
            Parameter            = p,
            Iso_type             = iso,
            Stratification       = "None",
            Population_Parameter = "mu.par",
            value                = vals[p]
          )
        }
      }
    }
    # Also add prec.logy draws (measurement precision)
    for (ch in seq_len(n_chains)) {
      for (it in seq_len(n_iter)) {
        # Precision ~ 5-15 on log scale
        pop_rows[[length(pop_rows) + 1]] <- tibble::tibble(
          Iteration            = it,
          Chain                = ch,
          Parameter            = "prec.logy",
          Iso_type             = iso,
          Stratification       = "None",
          Population_Parameter = "prec.logy",
          value                = runif(1, 5, 15)
        )
      }
    }
  }
  pop_df <- dplyr::bind_rows(pop_rows)

  # --- Construct the sr_model object ---
  structure(
    indiv_df,
    class             = c("sr_model", "tbl_df", "tbl", "data.frame"),
    nChains           = n_chains,
    nParameters       = as.integer(length(subjects) * 2 * 5),
    nIterations       = n_iter,
    nBurnin           = 100,
    nThin             = 1,
    population_params = pop_df,
    priors            = list(),
    fitted_residuals  = data.frame(),
    jags.post         = list()
  )
}


# ── Helper: create case_data object ───────────────────────────────────────

#' @param subjects_df Tibble with sid, cohort_name columns
#' @param antigen_key e.g. "IpaB", "Sf2a"
#' @param mfi_col_name e.g. "n_ipab_MFI"
make_mock_case_data <- function(subjects_df, antigen_key, mfi_col_name) {

  rows <- list()
  for (i in seq_len(nrow(subjects_df))) {
    sid <- subjects_df$sid[i]
    coh <- subjects_df$cohort_name[i]
    for (iso in isotypes) {
      for (tp_idx in seq_along(timepoints_nominal)) {
        tp <- timepoints_nominal[tp_idx]
        # Add small jitter to actual day
        actual_day <- tp + runif(1, -0.5, 0.5)
        if (actual_day < 1) actual_day <- tp
        # Simulate MFI: higher for matched serotype
        base_mfi <- if (grepl(tolower(substr(antigen_key, 1, 2)),
                              tolower(coh))) {
          exp(runif(1, 8, 12))  # matched: high signal
        } else {
          exp(runif(1, 3, 7))   # unmatched: low signal
        }
        # IpaB is always high (conserved)
        if (antigen_key == "IpaB") base_mfi <- exp(runif(1, 8, 12))
        # IgA decays faster
        if (iso == "IgA" && tp > 30) base_mfi <- base_mfi * 0.3

        rows[[length(rows) + 1]] <- tibble::tibble(
          isotype_name = iso,
          sid          = sid,
          timepoint    = tp,
          `Actual day` = round(actual_day, 2),
          cohort_name  = coh,
          !!mfi_col_name := round(base_mfi, 1),
          id           = sid,
          antigen_iso  = iso,
          visit        = tp,
          timeindays   = round(actual_day, 2),
          result       = round(base_mfi, 1),
          visit_num    = as.integer(tp_idx)
        )
      }
    }
  }

  out <- dplyr::bind_rows(rows)
  structure(
    out,
    class         = c("case_data", "tbl_df", "tbl", "data.frame"),
    id_var        = "id",
    biomarker_var = "antigen_iso",
    timeindays    = "timeindays",
    value_var     = "result"
  )
}


# ── Helper: create mock compiled data (raw Excel replacement) ─────────────

make_mock_compiled <- function() {
  rows <- list()
  for (i in seq_len(nrow(mock_subjects))) {
    sid <- mock_subjects$sid[i]
    coh <- mock_subjects$cohort_name[i]
    age <- mock_subjects$age[i]
    for (iso_num in 1:2) {
      iso_name <- isotypes[iso_num]
      for (tp in timepoints_nominal) {
        actual_day <- tp + runif(1, -0.5, 0.5)
        if (actual_day < 1) actual_day <- tp
        rows[[length(rows) + 1]] <- tibble::tibble(
          isotype       = iso_num,
          isotype_name  = iso_name,
          dilution      = 1,
          dilution_name = "1:1000",
          study_name    = "SOSAR",
          sampleID      = sid,
          sid           = sid,
          unq_id        = paste0(sid, "_d", tp, "_1:1000_", iso_name),
          cohort        = i,
          cohort_name   = coh,
          site          = 2,
          site_name     = "Dhaka",
          age           = age,
          treatment     = NA_character_,
          timepoint     = tp,
          `Actual day`  = round(actual_day, 2),
          ipab_MFI          = round(exp(runif(1, 8, 12)), 1),
          sf3aospbsa_MFI    = round(exp(runif(1, 4, 10)), 1),
          sf2aospbsa_MFI    = round(exp(runif(1, 4, 10)), 1),
          sf6ospbsa_MFI     = round(exp(runif(1, 3, 8)), 1),
          sonneiospbsa_MFI  = round(exp(runif(1, 4, 9)), 1),
          n_ipab_MFI        = round(exp(runif(1, 8, 12)), 1),
          n_sf3aospbsa_MFI  = round(exp(runif(1, 4, 10)), 1),
          n_sf2aospbsa_MFI  = round(exp(runif(1, 4, 10)), 1),
          n_sf6ospbsa_MFI   = round(exp(runif(1, 3, 8)), 1),
          n_sonneiospbsa_MFI = round(exp(runif(1, 4, 9)), 1),
          ipab_BC           = round(runif(1, 50, 200)),
          sf3aospbsa_BC     = round(runif(1, 50, 200)),
          sf2aospbsa_BC     = round(runif(1, 50, 200)),
          sf6ospbsa_BC      = round(runif(1, 50, 100)),
          sonneiospbsa_BC   = round(runif(1, 30, 80))
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}


# ══════════════════════════════════════════════════════════════════════════
# Generate all mock data objects
# ══════════════════════════════════════════════════════════════════════════

message("Generating mock data objects...")

all_sids   <- mock_subjects$sid
sf2a_sids  <- mock_subjects$sid[mock_subjects$cohort_name == "Sf2a"]
sonnei_sids <- mock_subjects$sid[mock_subjects$cohort_name == "sonnei"]
sf3a_sids  <- mock_subjects$sid[mock_subjects$cohort_name == "Sf3a"]
flexneri_sids <- mock_subjects$sid[mock_subjects$cohort_name %in%
                                     c("Sf2a", "Sf3a")]
under5_sids <- mock_subjects$sid[mock_subjects$age < 5]
plus5_sids  <- mock_subjects$sid[mock_subjects$age >= 5]

# --- Overall models (n = 5 mock subjects) ---
overall_IpaB_pop_6   <- make_mock_model(all_sids, "IpaB")
overall_Sf2a_pop_6   <- make_mock_model(all_sids, "Sf2a")
overall_Sf3a_pop_6   <- make_mock_model(all_sids, "Sf3a")
overall_Sf6_pop_6    <- make_mock_model(all_sids, "Sf6")
overall_Sonnei_pop_6 <- make_mock_model(all_sids, "Sonnei")

# --- Serotype-specific models ---
serotype_sf2a_3   <- make_mock_model(sf2a_sids, "Sf2a")
serotype_sf3a_3   <- make_mock_model(sf3a_sids, "Sf3a")
serotype_sonnei_3 <- make_mock_model(sonnei_sids, "Sonnei")

# --- Combined flexneri models ---
combined_flexneri_sf2a <- make_mock_model(flexneri_sids, "Sf2a")
combined_flexneri_sf3a <- make_mock_model(flexneri_sids, "Sf3a")

# --- Age-stratified IpaB models ---
overall_IpaB_pop_under5 <- make_mock_model(under5_sids, "IpaB")
overall_IpaB_pop_plus5  <- make_mock_model(plus5_sids,  "IpaB")

# --- Case data: overall (all 5 subjects) ---
dL_clean_Ipab_new   <- make_mock_case_data(mock_subjects, "IpaB",   "n_ipab_MFI")
dL_clean_sf2a_new   <- make_mock_case_data(mock_subjects, "Sf2a",   "n_sf2aospbsa_MFI")
dL_clean_sf3a_new   <- make_mock_case_data(mock_subjects, "Sf3a",   "n_sf3aospbsa_MFI")
dL_clean_sonnei_new <- make_mock_case_data(mock_subjects, "Sonnei", "n_sonneiospbsa_MFI")
dL_clean_sf6_new    <- make_mock_case_data(mock_subjects, "Sf6",    "n_sf6ospbsa_MFI")

# --- Case data: serotype-specific ---
dL_serotype_sf2a   <- make_mock_case_data(
  mock_subjects[mock_subjects$cohort_name == "Sf2a", ], "Sf2a", "n_sf2aospbsa_MFI"
)
dL_serotype_sf3a   <- make_mock_case_data(
  mock_subjects[mock_subjects$cohort_name == "Sf3a", ], "Sf3a", "n_sf3aospbsa_MFI"
)
dL_serotype_sonnei <- make_mock_case_data(
  mock_subjects[mock_subjects$cohort_name == "sonnei", ], "Sonnei", "n_sonneiospbsa_MFI"
)

# --- Case data: combined flexneri ---
dL_combined_flexneri_sf2a_case <- make_mock_case_data(
  mock_subjects[mock_subjects$cohort_name %in% c("Sf2a", "Sf3a"), ],
  "Sf2a", "n_sf2aospbsa_MFI"
)
dL_combined_flexneri_sf3a_case <- make_mock_case_data(
  mock_subjects[mock_subjects$cohort_name %in% c("Sf2a", "Sf3a"), ],
  "Sf3a", "n_sf3aospbsa_MFI"
)

# --- Case data: age-stratified ---
dL_clean_Ipab_new_under5 <- make_mock_case_data(
  mock_subjects[mock_subjects$age < 5, ], "IpaB", "n_ipab_MFI"
)
dL_clean_Ipab_new_plus5  <- make_mock_case_data(
  mock_subjects[mock_subjects$age >= 5, ], "IpaB", "n_ipab_MFI"
)

# --- Mock compiled data (replaces Excel file) ---
mock_compiled_data <- make_mock_compiled()

# --- Mock metadata ---
mock_metadata <- tibble::tibble(
  sampleID  = mock_subjects$sid,
  Gender    = sample(c("Male", "Female"), 5, replace = TRUE),
  HhIncome  = round(runif(5, 5000, 50000)),
  DiaBlood  = rep(1, 5),
  DiaWatery = sample(1:2, 5, replace = TRUE),
  Fev       = sample(1:2, 5, replace = TRUE),
  StoolNo   = round(runif(5, 4, 20)),
  MUAC      = round(runif(5, 12, 26), 1),
  IVFluidRec = rep(0, 5),
  HosDur    = round(runif(5, 0, 168))
)

# --- Mock duration of diarrhea ---
mock_durdia <- tibble::tibble(
  CaseID       = mock_subjects$sid,
  DurDia_hours = round(runif(5, 3, 96))
)


# ══════════════════════════════════════════════════════════════════════════
# Save all objects to data/
# ══════════════════════════════════════════════════════════════════════════

message("Saving mock data to data/ ...")

# Model objects
usethis::use_data(overall_IpaB_pop_6,   overwrite = TRUE)
usethis::use_data(overall_Sf2a_pop_6,   overwrite = TRUE)
usethis::use_data(overall_Sf3a_pop_6,   overwrite = TRUE)
usethis::use_data(overall_Sf6_pop_6,    overwrite = TRUE)
usethis::use_data(overall_Sonnei_pop_6, overwrite = TRUE)

usethis::use_data(serotype_sf2a_3,   overwrite = TRUE)
usethis::use_data(serotype_sf3a_3,   overwrite = TRUE)
usethis::use_data(serotype_sonnei_3, overwrite = TRUE)

usethis::use_data(combined_flexneri_sf2a, overwrite = TRUE)
usethis::use_data(combined_flexneri_sf3a, overwrite = TRUE)

usethis::use_data(overall_IpaB_pop_under5, overwrite = TRUE)
usethis::use_data(overall_IpaB_pop_plus5,  overwrite = TRUE)

# Case data objects
usethis::use_data(dL_clean_Ipab_new,   overwrite = TRUE)
usethis::use_data(dL_clean_sf2a_new,   overwrite = TRUE)
usethis::use_data(dL_clean_sf3a_new,   overwrite = TRUE)
usethis::use_data(dL_clean_sonnei_new, overwrite = TRUE)
usethis::use_data(dL_clean_sf6_new,    overwrite = TRUE)

usethis::use_data(dL_serotype_sf2a,   overwrite = TRUE)
usethis::use_data(dL_serotype_sf3a,   overwrite = TRUE)
usethis::use_data(dL_serotype_sonnei, overwrite = TRUE)

usethis::use_data(dL_combined_flexneri_sf2a_case, overwrite = TRUE)
usethis::use_data(dL_combined_flexneri_sf3a_case, overwrite = TRUE)

usethis::use_data(dL_clean_Ipab_new_under5, overwrite = TRUE)
usethis::use_data(dL_clean_Ipab_new_plus5,  overwrite = TRUE)

# Raw data replacements
usethis::use_data(mock_compiled_data, overwrite = TRUE)
usethis::use_data(mock_metadata,      overwrite = TRUE)
usethis::use_data(mock_durdia,        overwrite = TRUE)

message("Done! ", length(list.files("data", pattern = "\\.rda$")),
        " .rda files created in data/")
