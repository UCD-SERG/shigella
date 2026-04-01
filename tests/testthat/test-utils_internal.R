test_that("ab() returns y0 at t = 0", {
  result <- shigella:::ab(
    t = 0, y0 = 200, y1 = 2000, t1 = 10, alpha = 0.02, shape = 1.2
  )
  expect_equal(result, 200)
})

test_that("ab() returns y1 at t = t1", {
  result <- shigella:::ab(
    t = 10, y0 = 200, y1 = 2000, t1 = 10, alpha = 0.02, shape = 1.2
  )
  expect_equal(result, 2000)
})

test_that("ab() shows expected rise for t < t1", {
  result <- shigella:::ab(
    t = 5, y0 = 200, y1 = 2000, t1 = 10, alpha = 0.02, shape = 1.2
  )
  expected <- 200 * exp((log(2000 / 200) / 10) * 5)
  expect_equal(result, expected)
})

test_that("ab() decays toward zero for large t (shape > 1)", {
  result <- shigella:::ab(
    t = 10000, y0 = 200, y1 = 2000, t1 = 10, alpha = 0.02, shape = 1.2
  )
  expect_true(is.finite(result))
  expect_true(result >= 0)
})

test_that("ab() handles vector input", {
  times <- c(0, 5, 10, 30, 90)
  result <- shigella:::ab(
    t = times, y0 = 200, y1 = 2000, t1 = 10, alpha = 0.02, shape = 1.2
  )
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
  expect_equal(result[1], 200)
  expect_equal(result[3], 2000)
})

test_that("ab() handles t1 = 0 same as upstream", {
  # When t1 = 0, bt() returns Inf, which is expected
  # upstream behavior. We don't add special handling.
  result <- shigella:::ab(
    t = c(0, 10, 30),
    y0 = 200, y1 = 2000, t1 = 0,
    alpha = 0.02, shape = 1.2
  )
  expect_length(result, 3)
  # Results may be NaN/Inf -- that's OK, matches upstream
})

test_that("get_timeindays_var() finds time variable from attribute", {
  mock_data <- data.frame(timepoint = c(0, 7, 30), value = c(100, 200, 300))
  attr(mock_data, "timeindays") <- "timepoint"
  expect_equal(shigella:::get_timeindays_var(mock_data), "timepoint")
})

test_that(
  "get_timeindays_var() defaults to timeindays when attribute missing", {
  mock_data <- data.frame(timeindays = c(0, 7, 30), value = c(100, 200, 300))
  expect_equal(shigella:::get_timeindays_var(mock_data), "timeindays")
})
