test_that("write_status appends a line to the log file", {
  log_file <- tempfile(fileext = ".txt")
  on.exit(unlink(log_file))

  write_status(log_file, "TEST", "hello")

  lines <- readLines(log_file)
  expect_length(lines, 1L)
  expect_match(lines[1], "STEP=TEST")
  expect_match(lines[1], "hello")
})

test_that("write_status defaults msg to empty string", {
  log_file <- tempfile(fileext = ".txt")
  on.exit(unlink(log_file))

  write_status(log_file, "STEP_A")

  lines <- readLines(log_file)
  expect_length(lines, 1L)
  expect_match(lines[1], "STEP=STEP_A")
})

test_that("write_status appends multiple calls in order", {
  log_file <- tempfile(fileext = ".txt")
  on.exit(unlink(log_file))

  write_status(log_file, "FIRST")
  write_status(log_file, "SECOND")
  write_status(log_file, "THIRD")

  lines <- readLines(log_file)
  expect_length(lines, 3L)
  expect_match(lines[1], "FIRST")
  expect_match(lines[2], "SECOND")
  expect_match(lines[3], "THIRD")
})
