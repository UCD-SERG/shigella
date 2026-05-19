#' @keywords internal
#' @noRd
write_status <- function(status_file, step, msg = "") {
  cat(sprintf("[%s] STEP=%s | %s\n", format(Sys.time()), step, msg),
      file = status_file, append = TRUE)
}
