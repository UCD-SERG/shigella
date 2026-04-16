#' shigella: Antibody Kinetics and Serodynamics Tools for Shigella Research
#'
#' @description
#' The \pkg{shigella} package provides a complete, reproducible workflow for
#' the UCD-SERG Chapter 1 analysis of longitudinal *Shigella* antibody data.
#'
#' ## Core functions
#'
#' **Data processing**
#' - [process_shigella_data()] — filter and reshape raw Luminex data
#' - [load_all_antigens()] — convenience wrapper for all five antigens
#'
#' **Kinetic parameter extraction**
#' - [prep_pop_params()] — extract population-level MCMC draws
#' - [transform_pop_params()] — convert log-scale draws to natural scale
#' - [summarise_pop_params()] — posterior medians and 95% CrI
#' - [format_param_table()] — format for display
#' - [fmt_mci()] — format one median + CrI as a string
#'
#' **Posterior prediction**
#' - [predict_posterior_at_times()] — full posterior predictions
#' - [get_prediction_summary()] — summarised predictions (median + CrI)
#'
#' **Residual / accuracy metrics**
#' - [get_observed()] — extract observed data for one participant
#' - [compute_residual_metrics()] — MAE, RMSE, SSE at various summary levels
#' - [get_mae()] — per-individual MAE for one model-antigen-isotype
#'
#' **Tables**
#' - [build_table4()] — MAE comparison flextable across modelling approaches
#' - [build_kinetic_flextable()] — kinetic parameter flextable
#'
#' **Figures**
#' - [build_figure4_row()] — individual-level trajectory comparison
#' - [fig5_with_individuals()] — population trajectory with individual overlays
#' - [fig5_ipab_with_age()] — IpaB trajectories with age-stratified curves
#'
#' **Diagnostics**
#' - [posterior_pred()] — posterior predictive check plot
#' - [plot_model_comparison_forest()] — forest plot of kinetic parameters
#'
#' **Sensitivity analysis**
#' - [build_sensitivity_results()] — compute MAE across prior configurations
#' - [build_table_s3_gt()] — formatted sensitivity analysis table
#'
#' **Supplemental figures**
#' - [plot_age_spaghetti()] — age-stratified raw antibody trajectories
#' - [plot_mae_slopegraph()] — per-individual MAE slopegraph
#'
#' **Supplemental tables**
#' - [build_table_s4_gt()] — all-models parameter comparison
#' - [build_table_s5()] — raw per-individual MAE table
#'
#' ## Internal utilities
#' The package includes internal implementations of the two-phase antibody
#' kinetic model (`.ab()`, `.bt()`) and case-data helpers, eliminating
#' runtime dependencies on \pkg{serodynamics} and \pkg{serocalculator}
#' internals.
#'
#' ## Datasets
#' Pre-fitted model objects and case-data objects are stored in `data/`.
#' See [mock_compiled_data] and [overall_IpaB_pop_6] for documentation.
#' Mock versions are provided for testing; see `data-raw/mock_data.R`.
#'
#' @keywords internal
"_PACKAGE"
