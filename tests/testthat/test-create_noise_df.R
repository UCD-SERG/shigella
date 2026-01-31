test_that("create_noise_df creates correct structure", {
  result <- create_noise_df("MA USA")
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1)
})

test_that("create_noise_df has correct columns", {
  result <- create_noise_df("Ghana")
  
  expected_cols <- c("antigen_iso", "Country", "y.low", "eps", "nu", "y.high")
  expect_equal(names(result), expected_cols)
})

test_that("create_noise_df uses input country name", {
  countries <- c("MA USA", "Ghana", "Niger", "Sierra Leone")
  
  for (country in countries) {
    result <- create_noise_df(country)
    expect_equal(as.character(result$Country), country)
  }
})

test_that("create_noise_df has correct default values", {
  result <- create_noise_df("Test Country")
  
  expect_equal(result$antigen_iso, "IgG")
  expect_equal(result$y.low, 25)
  expect_equal(result$eps, 0.25)
  expect_equal(result$nu, 0.5)
  expect_equal(result$y.high, 200000)
})

test_that("create_noise_df Country is a factor", {
  result <- create_noise_df("Niger")
  
  expect_true(is.factor(result$Country))
})
