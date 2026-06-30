## =============================================================================
##  data-raw/02_serotype_models.R
##  Serotype-specific models (matched infection only):
##    Sf2a (n = 17), Sf3a (n = 8), S. sonnei (n = 11).
##  Source manuscript: manuscript2.qmd  (seed 6, 10,000,000 iterations).
##
##  Inputs  (00_build_case_data.R): dL_serotype_{sf2a,sf3a,sonnei}
##    -> the NON-"_new" serotype objects, i.e. ORIGINAL timeindays + cohort
##       filter. This matches manuscript2.qmd, which loads dL_serotype_{ag}
##       (not the "_new" variant) to fit these models. (The DurDia-adjusted
##       dL_serotype_{ag}_new objects exist for the supplement, not for fitting.)
##  Outputs (manuscript_data_dir):  serotype_{sf2a,sf3a,sonnei}_3.rda
## =============================================================================

source("data-raw/_config.R")

devtools::load_all() # loads shigella R/ functions (fit_and_save, load_inputs, ...)
library(serodynamics)

load_inputs(c(
  "dL_serotype_sf2a",
  "dL_serotype_sf3a",
  "dL_serotype_sonnei"
), dir = manuscript_data_dir)

set.seed(6) # preserved from manuscript2.qmd (see reproducibility note in _config.R)

start_time <- Sys.time()

## Order matches manuscript2.qmd: sf2a -> sf3a -> sonnei
fit_and_save(dL_serotype_sf2a, "serotype_sf2a_3", settings = mcmc_main, priors = prior_settings, dir = manuscript_data_dir)
fit_and_save(dL_serotype_sf3a, "serotype_sf3a_3", settings = mcmc_main, priors = prior_settings, dir = manuscript_data_dir)
fit_and_save(dL_serotype_sonnei, "serotype_sonnei_3", settings = mcmc_main, priors = prior_settings, dir = manuscript_data_dir)

message(
  "02_serotype_models.R runtime: ",
  format(round(Sys.time() - start_time, 2))
)
