#' Write a step-status line to a diagnostic log file
#'
#' @param status_file path to the status log file (character)
#' @param step character label for the current step
#' @param msg optional message string (default `""`)
#' @keywords internal
#' @export
write_status <- function(status_file, step, msg = "") {
  cat(sprintf("[%s] STEP=%s | %s\n", format(Sys.time()), step, msg),
      file = status_file, append = TRUE)
}
