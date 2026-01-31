test_that("process_shigella_data filters by study correctly", {
  # Create mock data with multiple studies
  mock_data <- data.frame(
    study_name = c(rep("SOSAR", 6), rep("GEMS", 4)),
    isotype_name = rep("IgG", 10),
    sid = rep(1:5, 2),
    timepoint = rep(c("V1", "V2"), 5),
    `Actual day` = rep(c(0, 14), 5),
    ipab_MFI = rnorm(10, 1000, 200),
    check.names = FALSE
  )
  
  result <- process_shigella_data(
    data = mock_data,
    study_filter = "SOSAR",
    antigen = ipab_MFI
  )
  
  # Should only include SOSAR rows
  expect_equal(nrow(result), 6)
})

test_that("process_shigella_data creates correct output columns", {
  mock_data <- data.frame(
    study_name = rep("SOSAR", 4),
    isotype_name = rep(c("IgG", "IgA"), each = 2),
    sid = c(1, 1, 2, 2),
    timepoint = rep(c("V1", "V2"), 2),
    `Actual day` = c(0, 14, 0, 28),
    sf2a_MFI = rnorm(4, 800, 100),
    check.names = FALSE
  )
  
  result <- process_shigella_data(
    data = mock_data,
    study_filter = "SOSAR",
    antigen = sf2a_MFI
  )
  
  # Check required output columns exist
  expected_cols <- c("index_id", "antigen_iso", "visit", "timeindays", 
                     "result", "visit_num")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("process_shigella_data removes NA timeindays", {
  mock_data <- data.frame(
    study_name = rep("SOSAR", 6),
    isotype_name = rep("IgG", 6),
    sid = rep(1:3, 2),
    timepoint = rep(c("V1", "V2"), 3),
    `Actual day` = c(0, 14, NA, 28, 0, NA),
    ipab_MFI = rnorm(6, 1000, 200),
    check.names = FALSE
  )
  
  result <- process_shigella_data(
    data = mock_data,
    study_filter = "SOSAR",
    antigen = ipab_MFI
  )
  
  # Should only have 4 rows (2 NAs removed)
  expect_equal(nrow(result), 4)
  expect_true(all(!is.na(result$timeindays)))
})

test_that("process_shigella_data assigns visit numbers correctly", {
  mock_data <- data.frame(
    study_name = rep("GEMS", 8),
    isotype_name = rep("IgG", 8),
    sid = rep(c(1, 1, 2, 2), 2),
    timepoint = rep(c("V1", "V2", "V1", "V2"), 2),
    `Actual day` = c(0, 14, 0, 14, 0, 28, 0, 21),
    sf3a_MFI = rnorm(8, 1200, 150),
    check.names = FALSE
  )
  
  result <- process_shigella_data(
    data = mock_data,
    study_filter = "GEMS",
    antigen = sf3a_MFI
  )
  
  # Check visit_num is assigned (should be 1, 2, etc. for each participant)
  expect_true("visit_num" %in% names(result))
  expect_true(all(result$visit_num > 0))
})

test_that("process_shigella_data handles different antigens", {
  mock_data <- data.frame(
    study_name = rep("SOSAR", 4),
    isotype_name = rep("IgA", 4),
    sid = c(1, 1, 2, 2),
    timepoint = rep(c("V1", "V2"), 2),
    `Actual day` = c(0, 14, 0, 14),
    antigen_col1 = rnorm(4, 1000, 200),
    antigen_col2 = rnorm(4, 800, 150),
    check.names = FALSE
  )
  
  # Test with first antigen
  result1 <- process_shigella_data(
    data = mock_data,
    study_filter = "SOSAR",
    antigen = antigen_col1
  )
  
  # Test with second antigen
  result2 <- process_shigella_data(
    data = mock_data,
    study_filter = "SOSAR",
    antigen = antigen_col2
  )
  
  expect_equal(nrow(result1), 4)
  expect_equal(nrow(result2), 4)
  expect_false(identical(result1$result, result2$result))
})
