## =============================================================================
##  data-raw/00_build_case_data.R
##  Raw Excel  ->  all dL_* case-data objects used by the model scripts (01-05)
##  and by the figure/table functions.
##
##  Run FIRST, after editing paths in data-raw/_config.R.
##  Requires devtools::load_all() (for the R/ data-prep functions) + serocalculator.
##
##  TIME-IN-DAYS CONVENTION (confirmed against the shipped objects + manuscripts)
##  --------------------------------------------------------------------------
##    "_new" in the name  -> DurDia-adjusted timeindays
##                           (timeindays += DurDia_hours/24; see
##                            adjust_timeindays_durdia()).
##    base / non-"_new"   -> ORIGINAL timeindays (= Actual day), no shift.
##    The TWO combined-flexneri objects are the explicit exception: they are
##    named "_case" (not "_new") and use ORIGINAL timeindays.
##
##  cohort_name (infecting serotype) and age live IN the compiled sheet
##  , so process_shigella_data()
##  keeps them and the subset_*() helpers filter on them directly -- no separate
##  metadata join.
##
##  WHAT EACH SECTION PRODUCES
##  --------------------------------------------------------------------------
##    2. dL_clean_{ag}            base, original times      (figure observed pts)
##       dL_clean_{ag}_new        DurDia-adjusted           (FIT overall models, 01)
##    3. dL_serotype_{ag}         original times + cohort filter
##                                                          (FIT serotype models, 02/05)
##       dL_serotype_{ag}_new     adjusted + cohort filter  (supplement loads)
##    4. dL_combined_flexneri_{ag}_case
##                                ORIGINAL times + (Sf2a+Sf3a) filter
##                                                          (FIT combined models, 04)
##    5. dL_clean_{ag}_new_{under5,plus5}
##                                adjusted + age filter     (FIT age models, 03)
## =============================================================================

source("data-raw/_config.R")

# devtools::load_all()                # makes R/data_prep.R functions available
library(readxl)
library(dplyr)

# -----------------------------------------------------------------------------
# 0. Read raw inputs
# -----------------------------------------------------------------------------
## Compiled MFI table: carries study_name, sid, timepoint, `Actual day`,
## isotype_name, cohort_name (J), age (M), and the antigen MFI columns.
compiled <- readxl::read_excel(raw_files$compiled, sheet = raw_sheets$compiled)

## --- DurDia (diarrhea duration prior to presentation) ---  
## Sheet "Duration of symptoms"; columns CaseID + DurDia_hours; 48 participants.
## Verified against dL_clean_sf2a_new: DurDia_days = DurDia_hours/24 and
## timeindays_new = (Actual day) + DurDia_days  (e.g. 3h -> +0.125 d).
durdia <- readxl::read_excel(raw_files$durdia, sheet = raw_sheets$durdia) |>
  dplyr::transmute(
    CaseID       = .data$CaseID,
    DurDia_hours = suppressWarnings(as.numeric(.data$DurDia_hours))
  )

# -----------------------------------------------------------------------------
# 1. Antigen map  (object stem  ->  MFI column in `compiled`)
# -----------------------------------------------------------------------------
antigen_map <- c(
  Ipab   = "n_ipab_MFI",
  sf2a   = "n_sf2aospbsa_MFI",
  sf3a   = "n_sf3aospbsa_MFI",
  sf6    = "n_sf6ospbsa_MFI",
  sonnei = "n_sonneiospbsa_MFI"
)

# helper: save an object under a chosen name into manuscript_data_dir
save_named <- function(obj, name) {
  assign(name, obj)
  save(list = name,
       file = file.path(manuscript_data_dir, paste0(name, ".rda")),
       compress = "xz",
       envir = environment())
  message("saved: ", name)
  invisible(obj)
}

# helper: base build = reshape one antigen, ORIGINAL timeindays (keeps cohort_name/age)
build_base <- function(ag) {
  process_shigella_data(
    data    = compiled,
    antigen = !!rlang::sym(antigen_map[[ag]])
  )
}

# helper: DurDia-adjusted ("_new") build = base + time shift
build_new <- function(ag) {
  build_base(ag) |>
    adjust_timeindays_durdia(durdia)
}

# -----------------------------------------------------------------------------
# 2. Overall case data: base (original times) + DurDia-adjusted ("_new")
#     base      -> dL_clean_{ag}        (original times; figure observed points)
#     adjusted  -> dL_clean_{ag}_new    (FIT the overall models, script 01)
#
#   Note: both now carry cohort_name + age (kept from compiled). The original
#   bare dL_clean_{ag} did not; this is a harmless 2-column superset that does
#   not affect any model fit (prep_data_new uses only id/antigen_iso/time/result)
#   or any figure.
# -----------------------------------------------------------------------------
for (ag in names(antigen_map)) {
  save_named(build_base(ag), paste0("dL_clean_", ag))
  save_named(build_new(ag),  paste0("dL_clean_", ag, "_new"))
}

# -----------------------------------------------------------------------------
# 3. Serotype-specific case data  (matched infection only, filter cohort_name)
#     dL_serotype_{ag}      ORIGINAL times + cohort filter
#                           -> input for serotype models (manuscript2/02,
#                              manuscript5/05 serotype_data)
#     dL_serotype_{ag}_new  ADJUSTED times + cohort filter  (supplement loads)
# -----------------------------------------------------------------------------
serotype_match <- c(sf2a = "Sf2a", sf3a = "Sf3a", sonnei = "sonnei")  # cohort_name value
for (ag in names(serotype_match)) {
  save_named(
    build_base(ag) |> subset_infecting_serotype(serotype_match[[ag]]),
    paste0("dL_serotype_", ag)
  )
  save_named(
    build_new(ag) |> subset_infecting_serotype(serotype_match[[ag]]),
    paste0("dL_serotype_", ag, "_new")
  )
}

# -----------------------------------------------------------------------------
# 4. Combined S. flexneri pool (Sf2a + Sf3a, n = 25), for the SF2a and SF3a beads
#     ORIGINAL timeindays (NO DurDia shift) -- confirmed: these "_case" objects
#     reproduce the manuscript4 fits, which use the unadjusted times.
# -----------------------------------------------------------------------------
for (ag in c("sf2a", "sf3a")) {
  save_named(
    build_base(ag) |> subset_combined_flexneri(),
    paste0("dL_combined_flexneri_", ag, "_case")
  )
}

# -----------------------------------------------------------------------------
# 5. Age-stratified case data (adjusted), all 5 antigens x {under5, plus5}
#     under5 = age < 5 ; plus5 = age > 5 ; age == 5 is in NEITHER group.
# -----------------------------------------------------------------------------
for (ag in names(antigen_map)) {
  new <- build_new(ag)
  save_named(subset_age_group(new, "under5"), paste0("dL_clean_", ag, "_new_under5"))
  save_named(subset_age_group(new, "plus5"),  paste0("dL_clean_", ag, "_new_plus5"))
}

message("\n00_build_case_data.R complete - all dL_* objects written to ",
        manuscript_data_dir)
