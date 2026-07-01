## =============================================================================
##  data-raw/04_combined_flexneri_models.R
##  Combined S. flexneri pool (Sf2a + Sf3a, n = 25), fit for the SF2a and SF3a
##  beads to exploit O-antigen cross-reactivity.
##  Source manuscript: manuscript4.qmd  (seed 9, 10,000,000 iterations).
##
##  Inputs  (00_build_case_data.R): dL_combined_flexneri_{sf2a,sf3a}_case
##  Outputs (manuscript_data_dir):  combined_flexneri_{sf2a,sf3a}.rda
## =============================================================================

source("data-raw/_config.R")

devtools::load_all() # loads shigella R/ functions (fit_and_save, load_inputs, ...)
library(serodynamics)

load_inputs(c(
  "dL_combined_flexneri_sf2a_case",
  "dL_combined_flexneri_sf3a_case"
), dir = manuscript_data_dir)

set.seed(9) # preserved from manuscript4.qmd

start_time <- Sys.time()

## Order matches manuscript4.qmd: sf2a -> sf3a
fit_and_save(dL_combined_flexneri_sf2a_case, "combined_flexneri_sf2a", settings = mcmc_main, priors = prior_settings, dir = manuscript_data_dir)
fit_and_save(dL_combined_flexneri_sf3a_case, "combined_flexneri_sf3a", settings = mcmc_main, priors = prior_settings, dir = manuscript_data_dir)

message(
  "04_combined_flexneri_models.R runtime: ",
  format(round(Sys.time() - start_time, 2))
)
