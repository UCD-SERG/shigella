## =============================================================================
##  data-raw/05_sensitivity_models.R
##  Prior-robustness sensitivity analysis (S3 Table).
##  Source manuscript: manuscript5.qmd  (seed 11, scaled to 3,000,000 iterations).
##
##  Four biomarkers x three prior configs x two model types = 24 fits.
##    biomarkers : Sonnei IgG, Sonnei IgA, Sf3a IgG, Sf3a IgA
##    priors     : primary, diffuse (prec/4), informative (prec*4)
##    model_type : overall  (dL_clean_{ag}_new),  serotype (dL_serotype_{ag})
##
##  Inputs  : dL_clean_{sonnei,sf3a}_new, dL_serotype_{sonnei,sf3a}
##  Outputs : sensitivity_{Biomarker}_{prior}_{model}.rda
##            -> each contains an object named `fit_obj`
##               (the S3-table loader in R/tables_supp.R expects this name).
## =============================================================================

source("data-raw/_config.R")

# devtools::load_all()
library(serodynamics)

load_inputs(c(
  "dL_clean_sonnei_new", "dL_serotype_sonnei",
  "dL_clean_sf3a_new",   "dL_serotype_sf3a"
))

set.seed(11)   # preserved from manuscript5.qmd

## Biomarker -> (antigen, isotype, overall data, serotype data). Order matches
## manuscript5.qmd so the single-seed run sequence is identical.
biomarker_map <- list(
  list(biomarker = "Sonnei IgG", overall = dL_clean_sonnei_new, serotype = dL_serotype_sonnei),
  list(biomarker = "Sonnei IgA", overall = dL_clean_sonnei_new, serotype = dL_serotype_sonnei),
  list(biomarker = "Sf3a IgG",   overall = dL_clean_sf3a_new,   serotype = dL_serotype_sf3a),
  list(biomarker = "Sf3a IgA",   overall = dL_clean_sf3a_new,   serotype = dL_serotype_sf3a)
)

prior_configs <- define_prior_configs()   # primary / diffuse / informative

start_time <- Sys.time()

for (bm in biomarker_map) {
  for (prior_name in names(prior_configs)) {
    for (model_type in c("overall", "serotype")) {

      dataset <- if (model_type == "overall") bm$overall else bm$serotype

      fit_name <- paste(gsub(" ", "_", bm$biomarker), prior_name, model_type, sep = "_")
      cat("Running:", fit_name, "\n")

      fit_and_save(
        data        = dataset,
        name        = paste0("sensitivity_", fit_name),
        object_name = "fit_obj",                 # S3-table loader contract
        settings    = mcmc_sensitivity,
        priors      = prior_configs[[prior_name]]
      )
    }
  }
}

message("05_sensitivity_models.R runtime: ",
        format(round(Sys.time() - start_time, 2)))
