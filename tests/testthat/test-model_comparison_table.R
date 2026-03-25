test_that("model_comparison_table computes differences correctly", {
  # Create mock metrics for two models
  metrics_overall <- data.frame(
    MAE = 150,
    RMSE = 200,
    SSE = 40000,
    n_obs = 100
  )
  
  metrics_pointwise <- data.frame(
    MAE = 120,
    RMSE = 180,
    SSE = 32400,
    n_obs = 100
  )
  
  result <- model_comparison_table(
    metrics_overall = metrics_overall,
    metrics_pointwise = metrics_pointwise,
    model_overall_label = "Overall Model",
    model_pointwise_label = "Pointwise Model"
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("Model", "MAE", "RMSE", "SSE", "n_obs") %in% names(result)))
  expect_true("delta_MAE" %in% names(result))
  expect_true("delta_RMSE" %in% names(result))
  
  # Check that differences are computed correctly
  delta_mae <- result$delta_MAE[result$Model == "Pointwise Model"]
  expect_equal(delta_mae, 120 - 150)
  
  delta_rmse <- result$delta_RMSE[result$Model == "Pointwise Model"]
  expect_equal(delta_rmse, 180 - 200)
})

test_that("model_comparison_table computes percent improvement", {
  metrics_overall <- data.frame(
    MAE = 100,
    RMSE = 200,
    SSE = 40000,
    n_obs = 100
  )
  
  metrics_pointwise <- data.frame(
    MAE = 80,
    RMSE = 160,
    SSE = 25600,
    n_obs = 100
  )
  
  result <- model_comparison_table(
    metrics_overall = metrics_overall,
    metrics_pointwise = metrics_pointwise
  )
  
  expect_true("pct_improve_MAE" %in% names(result))
  expect_true("pct_improve_RMSE" %in% names(result))
  
  # Check percent improvements are computed correctly
  # Pointwise model: (80 - 100) / 100 * 100 = -20% (20% improvement)
  pct_mae <- result$pct_improve_MAE[result$Model == "Pointwise Model"]
  expect_equal(pct_mae, -20)
  
  pct_rmse <- result$pct_improve_RMSE[result$Model == "Pointwise Model"]
  expect_equal(pct_rmse, -20)
})
