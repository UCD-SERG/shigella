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

test_that("prep_priors_stan includes biomarker priors only for model_2", {
  priors_model_2 <- prep_priors_stan(model = "model_2")
  priors_model_1 <- prep_priors_stan(model = "model_1")

  expect_true("tau_B_scale" %in% names(priors_model_2))
  expect_true("lkj_B_eta" %in% names(priors_model_2))
  expect_false("tau_B_scale" %in% names(priors_model_1))
  expect_false("lkj_B_eta" %in% names(priors_model_1))
})

test_that("prep_priors_stan model_2 structure is stable", {
  priors <- prep_priors_stan(model = "model_2")
  expect_equal(
    names(priors),
    c("mu_hyp_mean", "mu_hyp_sd", "tau_P_scale", "tau_eps_scale",
      "lkj_P_eta", "lkj_eps_eta", "tau_B_scale", "lkj_B_eta")
  )
  expect_true(all(vapply(priors, is.numeric, logical(1L))))
})
