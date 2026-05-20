#' Write a step-status line to a diagnostic log file
#'
#' Appends a timestamped status entry to a log file. Intended for use in
#' long-running HPC diagnostic scripts so crash location is visible even
#' when the script dies mid-way.
#'
#' @param status_file Path to the status log file (character).
#' @param step Character label for the current step.
#' @param msg Optional message string (default `""`).
#' @export
write_status <- function(status_file, step, msg = "") {
  cat(sprintf("[%s] STEP=%s | %s\n", format(Sys.time()), step, msg),
      file = status_file, append = TRUE)
}
