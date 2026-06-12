## =============================================================================
##  data-raw/03_age_models.R
##  Age-stratified models for every antigen x {under5, plus5}.
##  Source manuscript: manuscript3.qmd  (seed 8, 10,000,000 iterations).
##
##  Inputs  (00_build_case_data.R): dL_clean_{ag}_new_{under5,plus5}
##  Outputs (manuscript_data_dir):  overall_{Ag}_pop_{under5,plus5}.rda
##
##  NOTE: Only the IpaB age-stratified fits are used in the manuscript
##  (Table 3 / Fig 5A); the remaining antigen x age fits are produced here for
##  completeness because manuscript3.qmd fit all five under a single seed and the
##  run order affects nothing downstream once initsfunction() fixes the chains.
## =============================================================================

source("data-raw/_config.R")

# devtools::load_all()
library(serodynamics)

antigens <- c("Ipab", "sf2a", "sf3a", "sf6", "sonnei")
out_stem <- c(
  Ipab = "IpaB", sf2a = "Sf2a", sf3a = "Sf3a",
  sf6 = "Sf6", sonnei = "Sonnei"
) # output capitalization

load_inputs(as.vector(outer(
  paste0("dL_clean_", antigens, "_new"),
  c("_under5", "_plus5"),
  paste0
)))

set.seed(8) # preserved from manuscript3.qmd

start_time <- Sys.time()

## Order matches manuscript3.qmd: for each antigen, under5 then plus5.
for (ag in antigens) {
  fit_and_save(
    get(paste0("dL_clean_", ag, "_new_under5")),
    paste0("overall_", out_stem[[ag]], "_pop_under5"),
    settings = mcmc_main
  )
  fit_and_save(
    get(paste0("dL_clean_", ag, "_new_plus5")),
    paste0("overall_", out_stem[[ag]], "_pop_plus5"),
    settings = mcmc_main
  )
}

message(
  "03_age_models.R runtime: ",
  format(round(Sys.time() - start_time, 2))
)
