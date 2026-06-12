# Needed for data.table-style subset shig_bg[isotype == "IgG"] used in the plotting helper below.
utils::globalVariables("isotype")

#' Cross-reactivity heatmap (IgG | IgA, independent colour scales)
#'
#' @param compiled The `Compiled` sheet of the raw Shigella Excel.
#' @return A patchwork object.
#' @export
figure_crossreactivity_heatmap <- function(compiled) {
  shig_bg <- prep_heatmap_data(compiled)
  lut <- attr(shig_bg, "y_label_lookup")

  p_igg <- make_heatmap(shig_bg[isotype == "IgG"], "IgG", TRUE,  lut)
  p_iga <- make_heatmap(shig_bg[isotype == "IgA"], "IgA", FALSE, lut)

  p_igg + p_iga +
    patchwork::plot_layout(widths = c(1.1, 1)) +
    patchwork::plot_annotation(
      caption = paste("Each isotype panel uses an independent color scale.",
                      "Within each serotype facet, participants are sorted by age (oldest at top).", # nolint: line_length_linter.
                      sep = "\n"),
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(hjust = 0, size = 8, color = "grey40")) # nolint: line_length_linter.
    )
}
