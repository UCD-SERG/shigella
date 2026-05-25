test_that("compute_kinetics_at_time aborts on degenerate shape == 1", {
  expect_error(
    shigella:::.compute_kinetics_at_time(
      log_y0 = log(1), log_y1m0 = log(9), log_t1 = log(2),
      log_alpha = log(0.5), log_rm1 = -Inf, tt = 5
    ),
    regexp = "degenerate"
  )
})

test_that("compute_kinetics_at_time succeeds in growth phase even when log_rm1 = -Inf", {
  expect_no_error(
    shigella:::.compute_kinetics_at_time(
      log_y0 = log(1), log_y1m0 = log(9), log_t1 = log(2),
      log_alpha = log(0.5), log_rm1 = -Inf, tt = 1  # tt < t1_j = 2
    )
  )
})
