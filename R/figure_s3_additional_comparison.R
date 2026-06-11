#' S3 Fig: additional individual model-comparison rows
#'
#' Same three-row layout/colour coding as Figure 4, for different participants.
#' Wrapper over [figure_individual_comparison()].
#'
#' @param models,datasets Named lists as documented in
#'   [figure_individual_comparison()].
#' @param sids Length-3 subject ids (Sf2a, Sonnei, Sf3a); defaults are the
#'   supplement's choices.
#' @return A patchwork figure.
#' @export
figure_s3_additional_comparison <- function(models, datasets, # nolint: object_length_linter
                                            sids = c("SOSAR-23004", "SOSAR-21018",
                                                     "SOSAR-21002")) {
  figure_individual_comparison(models, datasets, sids = sids)
}
