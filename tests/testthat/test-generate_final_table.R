test_that("generate_final_table returns correct structure", {
  # Create mock simulation results
  mock_results <- lapply(1:200, function(i) {
    mock_est <- list()
    # Mock summary output
    mock_summary <- data.frame(
      incidence.rate = runif(1, 0.1, 0.3),
      SE = runif(1, 0.01, 0.05),
      CI.lwr = runif(1, 0.05, 0.15),
      CI.upr = runif(1, 0.2, 0.4)
    )
    
    # Create object with summary method
    mock_est$est1 <- structure(
      list(data = mock_summary),
      class = "mock_est"
    )
    mock_est
  })
  
  # Mock summary method for our test objects
  summary.mock_est <- function(object) {
    object$data
  }
  
  result <- generate_final_table(
    results_list = mock_results,
    sample_size = 200
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 200)
})

test_that("generate_final_table has correct columns", {
  # Create minimal mock results
  mock_results <- lapply(1:200, function(i) {
    list(
      est1 = structure(
        list(data = data.frame(
          incidence.rate = 0.15,
          SE = 0.03,
          CI.lwr = 0.09,
          CI.upr = 0.21
        )),
        class = "mock_est"
      )
    )
  })
  
  summary.mock_est <- function(object) {
    object$data
  }
  
  result <- generate_final_table(
    results_list = mock_results,
    sample_size = 150
  )
  
  expected_cols <- c("incidence.rate", "SE", "CI.lwr", "CI.upr", "index", "sample_size")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("generate_final_table includes sample size", {
  mock_results <- lapply(1:200, function(i) {
    list(
      est1 = structure(
        list(data = data.frame(
          incidence.rate = 0.15,
          SE = 0.03,
          CI.lwr = 0.09,
          CI.upr = 0.21
        )),
        class = "mock_est"
      )
    )
  })
  
  summary.mock_est <- function(object) {
    object$data
  }
  
  test_size <- 250
  result <- generate_final_table(
    results_list = mock_results,
    sample_size = test_size
  )
  
  expect_true(all(result$sample_size == test_size))
})

test_that("generate_final_table assigns index correctly", {
  mock_results <- lapply(1:200, function(i) {
    list(
      est1 = structure(
        list(data = data.frame(
          incidence.rate = runif(1, 0.1, 0.3),
          SE = 0.03,
          CI.lwr = 0.09,
          CI.upr = 0.21
        )),
        class = "mock_est"
      )
    )
  })
  
  summary.mock_est <- function(object) {
    object$data
  }
  
  result <- generate_final_table(
    results_list = mock_results,
    sample_size = 200
  )
  
  expect_equal(result$index, 1:200)
})
