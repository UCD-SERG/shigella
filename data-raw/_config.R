## =============================================================================
##  data-raw/_config.R
##  Shared settings for every model-generation script (01..05).
##  Source this at the top of each script:  source("data-raw/_config.R")
##
##  NOTHING here has side effects except defining objects + one helper.
## =============================================================================

## ====================== EDIT THESE AFTER DOWNLOADING DATA ====================
## (Ezra: point these at the raw files from SharePoint and a writable output dir.)
raw_data_dir        <- ""   # folder containing the raw Shigella Excel files
manuscript_data_dir <- ""   # output folder for fitted-model .rda files
## =============================================================================

stopifnot(
  "Set raw_data_dir in data-raw/_config.R"        = nzchar(raw_data_dir),
  "Set manuscript_data_dir in data-raw/_config.R" = nzchar(manuscript_data_dir)
)
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
  metadata = file.path(raw_data_dir,
                       "Additional metadata for Luminex sample set 10.28.2024.xlsx"),
  durdia   = file.path(raw_data_dir,
                       "Duration of diarrhea prior to presentation.xlsx")
)
raw_sheets <- list(
  compiled = "Compiled",              # confirm against the workbook
  metadata = "Metadata",              # confirm (Table 1 only)
  durdia   = "Duration of symptoms"   # [CONFIRMED] cols: CaseID + DurDia_hours (n=48)
)

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
define_prior_configs <- function(base = prior_settings) {
  list(
    primary     = utils::modifyList(base, list(label = "Primary")),
    diffuse     = utils::modifyList(base, list(label = "Diffuse",
                                               prec_hyp_param = base$prec_hyp_param / 4)),
    informative = utils::modifyList(base, list(label = "Informative",
                                               prec_hyp_param = base$prec_hyp_param * 4))
  )
}

## ---- Fit one model and save it ----------------------------------------------
## Reproduces the original `obj <- run_mod_pop(...); save(obj, file = ...)` idiom.
##
##  - `name`        : the .rda file stem (file is <name>.rda in `dir`).
##  - `object_name` : the variable name the fit is stored UNDER inside the .rda.
##                    Defaults to `name`, so scripts 01-04 produce
##                    overall_IpaB_pop_6.rda containing `overall_IpaB_pop_6`,
##                    matching the manuscript `load()` calls.
##                    Script 05 (sensitivity) passes object_name = "fit_obj"
##                    because the S3-table loader (tables_supp.R) expects every
##                    sensitivity_*.rda to contain an object literally called
##                    `fit_obj`.
##
##  Reproducibility note: serodynamics::initsfunction() assigns FIXED per-chain
##  RNG seeds (1-4) and named generators, so the JAGS draws are reproducible from
##  identical input data alone. The per-script set.seed() calls are preserved for
##  fidelity but do not themselves drive the JAGS output.
##
##  Requires `serodynamics` to be attached (see decision D1 in REPRODUCIBILITY_PLAN).
fit_and_save <- function(data,
                         name,
                         object_name = name,
                         settings    = mcmc_main,
                         priors      = prior_settings,
                         dir         = manuscript_data_dir,
                         with_post   = TRUE) {

  obj <- run_mod_pop(
    data     = data,
    file_mod = serodynamics_example("model.jags"),
    nchain   = settings$nchain,
    nadapt   = settings$nadapt,
    nburn    = settings$nburn,
    nmc      = settings$nmc,
    niter    = settings$niter,
    mu_hyp_param        = priors$mu_hyp_param,
    prec_hyp_param      = priors$prec_hyp_param,
    omega_param         = priors$omega_param,
    wishdf_param        = priors$wishdf_param,
    prec_logy_hyp_param = priors$prec_logy_hyp_param,
    with_post = with_post
  )

  assign(object_name, obj)
  save(list = object_name,
       file = file.path(dir, paste0(name, ".rda")),
       compress = "xz",
       envir = environment())
  message("saved: ", file.path(dir, paste0(name, ".rda")))
  invisible(obj)
}

## ---- Load a set of dL_* inputs by stem from manuscript_data_dir -------------
## Each <stem>.rda (written by 00_build_case_data.R) contains one object named
## <stem>; this loads them into the caller's environment.
load_inputs <- function(stems, dir = manuscript_data_dir,
                         envir = parent.frame()) {
  for (s in stems) {
    f <- file.path(dir, paste0(s, ".rda"))
    if (!file.exists(f)) {
      stop("Missing input: ", f, "\n  -> run data-raw/00_build_case_data.R first.")
    }
    load(f, envir = envir)
  }
  invisible(stems)
}
