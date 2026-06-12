## =============================================================================
##  data-raw/01_overall_models.R
##  Overall (pooled, n = 48) two-phase hierarchical models for all 5 antigens.
##  Source manuscript: manuscript.qmd  (seed 7, 10,000,000 iterations).
##
##  Inputs  (built by 00_build_case_data.R): dL_clean_{Ipab,sf2a,sf3a,sf6,sonnei}_new
##  Outputs (manuscript_data_dir):           overall_{IpaB,Sf2a,Sf3a,Sf6,Sonnei}_pop_6.rda
##
##  Run order is preserved exactly as in the manuscript so that, together with the
##  fixed per-chain RNG in serodynamics::initsfunction(), the fits are reproducible.
## =============================================================================

source("data-raw/_config.R")

# devtools::load_all()   # serodynamics + R/ helpers (run_mod_pop, fit_and_save, ...)
library(serodynamics)

## ---- Load DurDia-adjusted ("_new") overall case data ------------------------
load_inputs(c(
  "dL_clean_Ipab_new",
  "dL_clean_sf2a_new",
  "dL_clean_sonnei_new",
  "dL_clean_sf6_new",
  "dL_clean_sf3a_new"
))

## ---- Fit ---------------------------------------------------------------------
## set.seed() preserved from manuscript.qmd for fidelity (reproducibility is
## actually guaranteed by initsfunction()'s fixed per-chain seeds; see _config.R).
set.seed(7)

start_time <- Sys.time()

## Order matches manuscript.qmd: IpaB -> Sf2a -> Sonnei -> Sf6 -> Sf3a
fit_and_save(dL_clean_Ipab_new, "overall_IpaB_pop_6", settings = mcmc_main)
fit_and_save(dL_clean_sf2a_new, "overall_Sf2a_pop_6", settings = mcmc_main)
fit_and_save(dL_clean_sonnei_new, "overall_Sonnei_pop_6", settings = mcmc_main)
fit_and_save(dL_clean_sf6_new, "overall_Sf6_pop_6", settings = mcmc_main)
fit_and_save(dL_clean_sf3a_new, "overall_Sf3a_pop_6", settings = mcmc_main)

message(
  "01_overall_models.R runtime: ",
  format(round(Sys.time() - start_time, 2))
)
