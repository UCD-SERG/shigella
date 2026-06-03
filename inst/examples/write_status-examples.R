## Example: write_status()
##
## Appends a timestamped step-status line to a log file.

if (interactive()) {

  log_file <- tempfile(fileext = ".txt")

  write_status(log_file, "INIT",    "starting workflow")
  write_status(log_file, "STEP_1",  "data loaded")
  write_status(log_file, "STEP_2")  # msg defaults to ""
  write_status(log_file, "DONE",    "workflow complete")

  readLines(log_file)
  unlink(log_file)

}
