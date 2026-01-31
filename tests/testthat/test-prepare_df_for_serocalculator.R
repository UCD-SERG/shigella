test_that("prepare_df_for_serocalculator renames age column", {
  mock_df <- data.frame(
    id = 1:5,
    age_years = c(2, 5, 10, 15, 20),
    antibody = rnorm(5, 1000, 200)
  )
  
  result <- suppressMessages(
    prepare_df_for_serocalculator(
      df = mock_df,
      age_col = "age_years",
      value_col = "antibody"
    )
  )
  
  expect_true("age" %in% names(result))
  expect_false("age_years" %in% names(result))
})

test_that("prepare_df_for_serocalculator sets attributes", {
  mock_df <- data.frame(
    id = 1:5,
    age = c(2, 5, 10, 15, 20),
    value = rnorm(5, 1000, 200)
  )
  
  result <- suppressMessages(
    prepare_df_for_serocalculator(
      df = mock_df,
      age_col = "age",
      value_col = "value"
    )
  )
  
  expect_equal(attr(result, "age_var"), "age")
  expect_equal(attr(result, "value_var"), "value")
})

test_that("prepare_df_for_serocalculator handles default parameters", {
  mock_df <- data.frame(
    id = 1:5,
    age = c(2, 5, 10, 15, 20),
    value = rnorm(5, 1000, 200)
  )
  
  result <- suppressMessages(
    prepare_df_for_serocalculator(df = mock_df)
  )
  
  expect_true("age" %in% names(result))
  expect_equal(attr(result, "age_var"), "age")
  expect_equal(attr(result, "value_var"), "value")
})

test_that("prepare_df_for_serocalculator preserves other columns", {
  mock_df <- data.frame(
    id = 1:5,
    age = c(2, 5, 10, 15, 20),
    country = rep("Ghana", 5),
    value = rnorm(5, 1000, 200),
    study = rep("SOSAR", 5)
  )
  
  result <- suppressMessages(
    prepare_df_for_serocalculator(
      df = mock_df,
      age_col = "age",
      value_col = "value"
    )
  )
  
  expect_true(all(c("id", "country", "study") %in% names(result)))
  expect_equal(nrow(result), 5)
})
