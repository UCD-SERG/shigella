test_that("plot_summary_metrics returns ggplot object", {
  mock_metrics <- data.frame(
    sample_size = rep(c(100, 200, 300), each = 3),
    empirical_se = runif(9, 0.01, 0.1),
    Age_Group = rep(c("<5", "5-15", "16+"), 3)
  )
  
  result <- plot_summary_metrics(mock_metrics)
  
  expect_s3_class(result, "ggplot")
})

test_that("plot_summary_metrics handles different age groups", {
  mock_metrics <- data.frame(
    sample_size = rep(c(100, 200), each = 2),
    empirical_se = c(0.05, 0.08, 0.03, 0.06),
    Age_Group = rep(c("<5", "16+"), 2)
  )
  
  result <- plot_summary_metrics(mock_metrics)
  
  expect_s3_class(result, "ggplot")
  # Plot should be created without errors
  expect_true(TRUE)
})

test_that("plot_summary_metrics requires correct columns", {
  # Missing required column
  mock_metrics_bad <- data.frame(
    sample_size = c(100, 200, 300),
    wrong_col = runif(3, 0.01, 0.1),
    Age_Group = c("<5", "5-15", "16+")
  )
  
  # Should error when empirical_se is missing
  expect_error(plot_summary_metrics(mock_metrics_bad))
})
