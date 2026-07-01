## =============================================================================
##  data-raw/_config.R
##  Shared settings for every model-generation script (01..05).
##  Source this at the top of each script:  source("data-raw/_config.R")
##
##  This file defines settings objects only (paths, MCMC parameters, priors,
##  sheet names).  Function definitions (define_prior_configs, fit_and_save,
##  load_inputs) live in R/ and are loaded via devtools::load_all() or
##  library(shigella) in each script.
## =============================================================================

## ====================== EDIT THESE AFTER DOWNLOADING DATA ====================
## (Ezra: point these at the raw files from SharePoint and a writable output dir.)
raw_data_dir <- "" # folder containing the raw Shigella Excel files
manuscript_data_dir <- "" # output folder for fitted-model .rda files
## =============================================================================

if (!nzchar(raw_data_dir)) {
  cli::cli_abort("Set {.code raw_data_dir} in {.file data-raw/_config.R}")
}
if (!nzchar(manuscript_data_dir)) {
  cli::cli_abort("Set {.code manuscript_data_dir} in {.file data-raw/_config.R}")
}
if (!dir.exists(manuscript_data_dir)) {
  dir.create(manuscript_data_dir, recursive = TRUE)
}

## ---- Raw Excel filenames (confirm against SharePoint) -----------------------
## NOTE on what each workbook supplies:
##   * compiled : the MFI data + sid/timepoint/`Actual day`/isotype_name AND
##                cohort_name (infecting serotype) + age  -> drives 00_build and
##                every model. (cohort_name/age are columns IN this sheet.)
##   * metadata : Table 1 CLINICAL variables only (Gender, DiaBlood, DiaWatery,
##                Fev, MUAC, HosDur). Not needed by 00_build or the models.
##   * durdia   : DurDia_hours per CaseID -> the "_new" time shift.
raw_files <- list(
  compiled = file.path(raw_data_dir, "3.8.2024 Compiled Shigella datav2.xlsx"),
  metadata = file.path(
    raw_data_dir,
    "Additional metadata for Luminex sample set 10.28.2024.xlsx"
  ),
  durdia = file.path(
    raw_data_dir,
    "Duration of diarrhea prior to presentation.xlsx"
  )
)
raw_sheets <- list(
  compiled = "Compiled", # confirm against the workbook
  metadata = "Metadata", # confirm (Table 1 only)
  durdia   = "Duration of symptoms" # [CONFIRMED] cols: CaseID + DurDia_hours (n=48)
)

## ---- Study-wide constants ---------------------------------------------------
## Days from symptom onset to the day-0 (enrollment) blood draw.
symptom_onset_offset_days <- 2L

## ---- MCMC settings ----------------------------------------------------------
## Main models (manuscript .qmd 1-4): 10,000,000 iterations.
mcmc_main <- list(
  nchain = 4,
  nadapt = 25000,
  nburn  = 50000,
  nmc    = 15000,
  niter  = 10000000
)

## Sensitivity analysis (manuscript5.qmd): scaled down to 3,000,000 iterations.
mcmc_sensitivity <- list(
  nchain = 4,
  nadapt = 7500,
  nburn  = 15000,
  nmc    = 4500,
  niter  = 3000000
)

## ---- Prior settings (primary specification, S1 Table) -----------------------
prior_settings <- list(
  mu_hyp_param        = c(8, 10, 0.5, -4, -1),
  prec_hyp_param      = c(0.25, 1e-4, 0.5, 0.001, 0.25),
  omega_param         = c(1, 1, 1, 1, 1),
  wishdf_param        = 8,
  prec_logy_hyp_param = c(1, 1)
)

## Alternative prior configurations for the sensitivity analysis.
## diffuse = prec/4 (wider), informative = prec*4 (tighter).
## Function define_prior_configs() is defined in R/define_prior_configs.R.
