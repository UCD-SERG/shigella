test_that("predict_posterior_at_times returns expected structure", {
  # Create minimal mock posterior data
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
  
  times <- c(0, 30, 90)
  
  result <- shigella:::predict_posterior_at_times(
    model = mock_model,
    ids = "test_id",
    antigen_iso = "IgG",
    times = times
  )
  
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(all(c("id", "t", "res") %in% names(result)))
  expect_true(all(result$t %in% times))
  expect_true(all(is.finite(result$res)))
  expect_true(all(result$res > 0))
})

test_that("predict_posterior_at_times handles multiple subjects", {
  # Create mock data for multiple subjects
  mock_model <- rbind(
    data.frame(
      Subject = rep("id1", 15),
      Iso_type = rep("IgG", 15),
      Chain = rep(1, 15),
      Iteration = rep(1:3, 5),
      Parameter = rep(c("y0", "y1", "t1", "alpha", "shape"), each = 3),
      value = c(200, 210, 190, 2000, 2100, 1900, 10, 12, 9, 
                0.02, 0.025, 0.015, 1.0, 1.1, 0.9)
    ),
    data.frame(
      Subject = rep("id2", 15),
      Iso_type = rep("IgG", 15),
      Chain = rep(1, 15),
      Iteration = rep(1:3, 5),
      Parameter = rep(c("y0", "y1", "t1", "alpha", "shape"), each = 3),
      value = c(250, 260, 240, 2500, 2600, 2400, 12, 14, 11, 
                0.03, 0.035, 0.025, 1.2, 1.3, 1.1)
    )
  )
  
  result <- shigella:::predict_posterior_at_times(
    model = mock_model,
    ids = c("id1", "id2"),
    antigen_iso = "IgG",
    times = c(0, 30)
  )
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("id1", "id2") %in% result$id))
  expect_true(all(is.finite(result$res)))
})
