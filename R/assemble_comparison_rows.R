#' Assemble a 3-row individual model-comparison figure
#'
#' @param rows A list of patchwork rows from [build_figure4_row()].
#' @return A patchwork figure with a shared legend strip.
#' @export
assemble_comparison_rows <- function(rows) {
  body <- Reduce(`/`, rows)
  (body / patchwork::wrap_elements(model_comparison_legend())) +
    patchwork::plot_layout(heights = c(rep(1, length(rows)), 0.15))
}
