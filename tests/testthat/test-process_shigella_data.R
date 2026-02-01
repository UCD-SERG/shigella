test_that("process_shigella_data filters and reshapes correctly", {
  # Create mock input data
  mock_raw_data <- data.frame(
    study_name = rep(c("SOSAR", "OTHER"), each = 6),
    isotype_name = rep(c("IgG", "IgA"), 6),
    sid = rep(c("ID001", "ID002"), each = 6),
    timepoint = rep(c("V1", "V2", "V3"), 4),
    `Actual day` = rep(c(0, 7, 30), 4),
    n_ipab_MFI = runif(12, 100, 5000),
    check.names = FALSE
  )
  
  result <- process_shigella_data(
    data = mock_raw_data,
    study_filter = "SOSAR",
    antigen = n_ipab_MFI
  )
  
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(all(c("index_id", "antigen_iso", "visit", "timeindays", "result") %in% names(result)))
  expect_true(all(result$index_id %in% c("ID001", "ID002")))
  expect_true(all(result$antigen_iso %in% c("IgG", "IgA")))
  expect_true(all(!is.na(result$timeindays)))
})

test_that("process_shigella_data removes NA timepoints", {
  mock_raw_data <- data.frame(
    study_name = rep("SOSAR", 6),
    isotype_name = rep("IgG", 6),
    sid = rep("ID001", 6),
    timepoint = c("V1", "V2", "V3", "V4", "V5", "V6"),
    `Actual day` = c(0, 7, NA, 30, NA, 90),
    n_ipab_MFI = runif(6, 100, 5000),
    check.names = FALSE
  )
  
  result <- process_shigella_data(
    data = mock_raw_data,
    study_filter = "SOSAR",
    antigen = n_ipab_MFI
  )
  
  expect_true(all(!is.na(result$timeindays)))
  expect_equal(nrow(result), 4)  # Should only have 4 rows (excluding 2 NAs)
})
