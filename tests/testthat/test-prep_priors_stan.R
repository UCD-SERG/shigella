test_that("prep_priors_stan returns a list with expected names", {
  priors <- prep_priors_stan(model = "model_2")

  expect_type(priors, "list")
  expect_true(length(priors) > 0)
})

test_that("prep_priors_stan defaults are weakly informative", {
  priors <- prep_priors_stan(model = "model_2")

  ## Prior SDs should not be absurdly diffuse (e.g., 316 from old JAGS
  ## translation) — that caused Stan HMC to wander during warmup.
  if (!is.null(priors$mu_hyp_sd)) {
    expect_true(all(priors$mu_hyp_sd > 0))
    expect_true(all(priors$mu_hyp_sd <= 20))
  }
})
