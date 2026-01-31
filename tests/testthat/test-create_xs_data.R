test_that("create_xs_data filters data correctly", {
  # Create mock data
  mock_df <- data.frame(
    sid = 1:10,
    site_name = c(rep("MA USA", 5), rep("Ghana", 5)),
    study_name = rep("SOSAR", 10),
    age = c(3, 8, 12, 18, 45, 2, 7, 11, 20, 50),
    isotype_name = rep("IgG", 10),
    n_ipab_MFI = rnorm(10, 1000, 200)
  )
  
  # Test filtering by country
  result <- create_xs_data(
    df = mock_df,
    filter_countries = "MA USA",
    filter_antigen_iso = "IgG",
    value_col = n_ipab_MFI
  )
  
  expect_equal(nrow(result), 5)
  expect_true(all(result$Country == "MA USA"))
  expect_true(all(result$antigen_iso == "IgG"))
})

test_that("create_xs_data creates age categories correctly", {
  mock_df <- data.frame(
    sid = 1:6,
    site_name = rep("Ghana", 6),
    study_name = rep("SOSAR", 6),
    age = c(3, 8, 17, 25, 50, 1),
    isotype_name = rep("IgG", 6),
    n_ipab_MFI = rnorm(6, 1000, 200)
  )
  
  result <- create_xs_data(
    df = mock_df,
    filter_countries = "Ghana",
    filter_antigen_iso = "IgG",
    value_col = n_ipab_MFI
  )
  
  expect_true(all(c("<5", "5-15", "16+") %in% result$ageCat))
  expect_equal(sum(result$ageCat == "<5"), 2)
  expect_equal(sum(result$ageCat == "5-15"), 1)
  expect_equal(sum(result$ageCat == "16+"), 3)
})

test_that("create_xs_data handles multiple countries", {
  mock_df <- data.frame(
    sid = 1:15,
    site_name = rep(c("MA USA", "Ghana", "Niger"), each = 5),
    study_name = rep("SOSAR", 15),
    age = runif(15, 1, 60),
    isotype_name = rep("IgG", 15),
    n_ipab_MFI = rnorm(15, 1000, 200)
  )
  
  result <- create_xs_data(
    df = mock_df,
    filter_countries = c("MA USA", "Ghana"),
    filter_antigen_iso = "IgG",
    value_col = n_ipab_MFI
  )
  
  expect_equal(nrow(result), 10)
  expect_true(all(result$Country %in% c("MA USA", "Ghana")))
})

test_that("create_xs_data returns correct column names", {
  mock_df <- data.frame(
    sid = 1:5,
    site_name = rep("Niger", 5),
    study_name = rep("GEMS", 5),
    age = c(2, 5, 10, 15, 20),
    isotype_name = rep("IgA", 5),
    n_sf2a_MFI = rnorm(5, 800, 150)
  )
  
  result <- create_xs_data(
    df = mock_df,
    filter_countries = "Niger",
    filter_antigen_iso = "IgA",
    value_col = n_sf2a_MFI
  )
  
  expected_cols <- c("id", "Country", "study", "age", "antigen_iso", "value", "ageCat")
  expect_equal(names(result), expected_cols)
})
