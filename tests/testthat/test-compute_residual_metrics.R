test_that("compute_residual_metrics returns expected structure", {
  # Note: This test uses mock data structure
  # Once mock_posterior_draws and mock_case_data are generated,
  # these tests should work
  
  # Create minimal mock data inline for testing
  mock_model <- data.frame(
    Subject = rep("test_id", 15),
    Iso_type = rep("IgG", 15),
    Chain = rep(1, 15),
    Iteration = rep(1:3, 5),
    Parameter = rep(c("y0", "y1", "t1", "alpha", "shape"), each = 3),
    value = c(
      200, 210, 190,  # y0
      2000, 2100, 1900,  # y1
      10, 12, 9,  # t1
      0.02, 0.025, 0.015,  # alpha
      1.0, 1.1, 0.9  # shape
    )
  )
  
  mock_dataset <- data.frame(
    id = rep("test_id", 4),
    antigen_iso = rep("IgG", 4),
    timepoint = c(0, 7, 30, 90),
    value = c(250, 1500, 1200, 800)
  )
  attr(mock_dataset, "timeindays") <- "timepoint"
  attr(mock_dataset, "value_var") <- "value"
  class(mock_dataset) <- c("case_data", "data.frame")
  
  # Test pointwise summary
  result_pw <- compute_residual_metrics(
    model = mock_model,
    dataset = mock_dataset,
    ids = "test_id",
    antigen_iso = "IgG",
    scale = "original",
    summary_level = "pointwise"
  )
  
  expect_s3_class(result_pw, "data.frame")
  expect_true(nrow(result_pw) > 0)
  expect_true(all(c("id", "antigen_iso", "t", "obs", "pred_med", 
                    "residual", "abs_residual", "sq_residual") %in% names(result_pw)))
  expect_true(all(is.finite(result_pw$residual)))
  
  # Test id_antigen summary
  result_id <- compute_residual_metrics(
    model = mock_model,
    dataset = mock_dataset,
    ids = "test_id",
    antigen_iso = "IgG",
    scale = "original",
    summary_level = "id_antigen"
  )
  
  expect_s3_class(result_id, "data.frame")
  expect_equal(nrow(result_id), 1)
  expect_true(all(c("id", "antigen_iso", "MAE", "RMSE", "SSE", "n_obs") %in% names(result_id)))
  expect_true(result_id$MAE >= 0)
  expect_true(result_id$RMSE >= 0)
  expect_true(result_id$SSE >= 0)
  expect_true(result_id$n_obs > 0)
})

test_that("compute_residual_metrics handles log scale", {
  # Minimal mock data
  mock_model <- data.frame(
    Subject = rep("test_id", 15),
    Iso_type = rep("IgG", 15),
    Chain = rep(1, 15),
    Iteration = rep(1:3, 5),
    Parameter = rep(c("y0", "y1", "t1", "alpha", "shape"), each = 3),
    value = c(200, 210, 190, 2000, 2100, 1900, 10, 12, 9, 
              0.02, 0.025, 0.015, 1.0, 1.1, 0.9)
  )
  
  mock_dataset <- data.frame(
    id = rep("test_id", 3),
    antigen_iso = rep("IgG", 3),
    timepoint = c(7, 30, 90),
    value = c(1500, 1200, 800)
  )
  attr(mock_dataset, "timeindays") <- "timepoint"
  attr(mock_dataset, "value_var") <- "value"
  class(mock_dataset) <- c("case_data", "data.frame")
  
  result_log <- compute_residual_metrics(
    model = mock_model,
    dataset = mock_dataset,
    ids = "test_id",
    antigen_iso = "IgG",
    scale = "log",
    summary_level = "overall"
  )
  
  expect_s3_class(result_log, "data.frame")
  expect_equal(nrow(result_log), 1)
  expect_true(all(is.finite(result_log$MAE)))
  expect_true(all(is.finite(result_log$RMSE)))
})
