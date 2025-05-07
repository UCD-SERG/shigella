test_that(
  desc = "process_shigella_data() returns case data structure without errors",
  code = {
    # Simulate input data directly instead of reading from xlsx file
    df <- tibble::tibble(
      study_name   = "SOSAR",
      isotype_name = c("IgA", "IgG"),
      sid          = c("ID1", "ID2"),
      timepoint    = c("D0", "D28"),
      `Actual day` = c(0, 28),
      n_ipab_MFI   = c(1000, 800)
    )
    
    # Suppress warnings and run the function
    result <- suppressWarnings(
      process_shigella_data(
        data = df,
        study_filter = "SOSAR",
        antigen = n_ipab_MFI
      )
    )
    
    # Assertions
    expect_s3_class(result, "case_data")
    expect_true(all(c("index_id", "antigen_iso", "timeindays", "result") %in% names(result)))
    
    # Optional: Snapshot test
    testthat::expect_snapshot(head(result))
  }
)

