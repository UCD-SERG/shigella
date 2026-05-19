# Helper: resolve and create the compile output directory.
# Priority: argument > environment variable > /tmp fallback.
#' @keywords internal
#' @noRd
.setup_compile_dir <- function(compile_dir) {
  if (is.null(compile_dir)) {
    compile_dir <- Sys.getenv("STAN_COMPILE_DIR", unset = "")
    if (compile_dir == "") {
      user <- Sys.getenv("USER", unset = "default")
      compile_dir <- file.path("/tmp", user, "cmdstan_bin")
    }
  }
  if (!dir.exists(compile_dir)) {
    dir.create(compile_dir, recursive = TRUE, mode = "0755")
  }
  compile_dir
}
