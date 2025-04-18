test_that(
  desc = "process_shigella_data() returns case data structure without errors",
  code = {
    
    # Load an example dataset from the package
    file_path <- system.file("extdata", "3.8.2024 Compiled Shigella datav2.xlsx", package = "shigella")
    
    # Read the Excel sheet
    df <- readxl::read_excel(file_path, sheet = "Compiled")
    
    # Suppress warnings (e.g., due to missing data)
    suppressWarnings({
      # Run the function
      result <- process_shigella_data(
        data = df,
        study_filter = "SOSAR",
        antigen = n_ipab_MFI
      )
    }) |> expect_no_error()
    
    # Ensure the result has expected structure (is case data and has known columns)
    expect_s3_class(result, "case_data")
    expect_true(all(c("index_id", "antigen_iso", "timeindays", "result") %in% names(result)))
    
    # Snapshot test to ensure consistency in structure and preview data
    testthat::expect_snapshot(head(result))
  }
)
