// ============================================================
// model_1.stan
//
// Model 1 (formerly Model A): Independent biomarker antibody kinetics.
// Each biomarker fit independently — no cross-biomarker correlation.
//
// Log-space kinetics: uses log_two_phase_curve() ported from model_2.stan.
// This removes the discontinuous fallback (no `if base <= 0`) and aligns
// the likelihood implementation with the validated model_2.stan.
//
// Compatible with serodynamics::prep_data_stan() output:
//   N, K, P, max_obs, n_obs[N], time_obs[N, max_obs], log_y[N, max_obs, K]
//
// Compatible with prep_priors_stan(model = "model_1") output.
// ============================================================

functions {
  // Compute log(y(t)) DIRECTLY (matching JAGS reference and model_2.stan)
  // Uses log-space parameters where natural:
  //   - log_y0, log_y1 are log of baseline / peak antibody
  //   - t1, alpha, shape are in natural scale
  real log_two_phase_curve(real t,
                            real log_y0, real log_y1,
                            real t1, real alpha, real shape) {
    if (t <= t1) {
      // Active phase: log(y(t)) = log(y0) + beta * t
      //   where beta = (log(y1) - log(y0)) / t1
      real beta = (log_y1 - log_y0) / t1;
      return log_y0 + beta * t;
    } else {
      // Recovery phase: log(y(t)) = 1/(1-shape) * log(inside)
      // inside = y1^(1-shape) - (1-shape)*alpha*(t-t1)
      //
      // Since shape > 1, (1 - shape) < 0:
      //   y1^(1-shape) = exp((1-shape) * log_y1)  (small positive)
      //   -(1-shape)*alpha*(t-t1) = (shape-1)*alpha*(t-t1)  (positive)
      // So inside = small_positive + positive = positive  ✓
      // No need for fallback because both terms are guaranteed positive.
      real one_minus_shape = 1 - shape;
      real first_term  = exp(one_minus_shape * log_y1);
      real second_term = (shape - 1) * alpha * (t - t1);
      real inside = first_term + second_term;
      return log(inside) / one_minus_shape;
    }
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;
  int<lower=1> max_obs;
  array[N] int<lower=0, upper=max_obs> n_obs;
  array[N, max_obs] real time_obs;
  array[N, max_obs, K] real log_y;

  vector[P] mu_hyp_mean;
  vector<lower=0>[P] mu_hyp_sd;
  real<lower=0> tau_P_scale;
  real<lower=0> tau_eps_scale;
  real<lower=0> lkj_P_eta;
  // lkj_eps_eta is NOT used here: model_1 uses independent per-biomarker
  // residual scales (tau_eps[k]), not an LKJ prior on epsilon correlation.
  // For the LKJ epsilon prior, see model_2.stan.
}

parameters {
  matrix[K, P] M;
  array[K] cholesky_factor_corr[P] L_Omega_P;
  array[K] vector<lower=0>[P] tau_P;
  vector<lower=0>[K] tau_eps;
  array[K] matrix[N, P] Z;
}

transformed parameters {
  array[N, K] vector[P] theta;
  for (k in 1:K) {
    matrix[P, P] L_Sigma_P_k = diag_pre_multiply(tau_P[k], L_Omega_P[k]);
    for (i in 1:N) {
      theta[i, k] = M[k]' + L_Sigma_P_k * Z[k][i]';
    }
  }
}

model {
  for (k in 1:K) {
    for (p in 1:P) {
      M[k, p] ~ normal(mu_hyp_mean[p], mu_hyp_sd[p]);
    }
    L_Omega_P[k] ~ lkj_corr_cholesky(lkj_P_eta);
    tau_P[k] ~ cauchy(0, tau_P_scale);
    to_vector(Z[k]) ~ std_normal();
  }
  tau_eps ~ cauchy(0, tau_eps_scale);

  for (i in 1:N) {
    if (n_obs[i] > 0) {
      for (o in 1:n_obs[i]) {
        for (k in 1:K) {
          real log_y0_o = theta[i, k][1];
          real log_y1_o = log_sum_exp(log_y0_o, theta[i, k][2]);
          real t1_o     = exp(theta[i, k][3]);
          real alpha_o  = exp(theta[i, k][4]);
          real shape_o  = exp(theta[i, k][5]) + 1;
          real mu_log = log_two_phase_curve(time_obs[i, o], log_y0_o, log_y1_o,
                                             t1_o, alpha_o, shape_o);
          log_y[i, o, k] ~ normal(mu_log, tau_eps[k]);
        }
      }
    }
  }
}

generated quantities {
  array[K] corr_matrix[P] Omega_P;
  for (k in 1:K) {
    Omega_P[k] = multiply_lower_tri_self_transpose(L_Omega_P[k]);
  }

  array[N, K] real y0;
  array[N, K] real y1;
  array[N, K] real t1;
  array[N, K] real alpha;
  array[N, K] real shape;
  for (i in 1:N) {
    for (k in 1:K) {
      y0[i, k]    = exp(theta[i, k][1]);
      y1[i, k]    = y0[i, k] + exp(theta[i, k][2]);
      t1[i, k]    = exp(theta[i, k][3]);
      alpha[i, k] = exp(theta[i, k][4]);
      shape[i, k] = exp(theta[i, k][5]) + 1;
    }
  }

  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = 0;
    if (n_obs[i] > 0) {
      for (o in 1:n_obs[i]) {
        for (k in 1:K) {
          real log_y0_o = theta[i, k][1];
          real log_y1_o = log_sum_exp(log_y0_o, theta[i, k][2]);
          real t1_o     = exp(theta[i, k][3]);
          real alpha_o  = exp(theta[i, k][4]);
          real shape_o  = exp(theta[i, k][5]) + 1;
          real mu_log = log_two_phase_curve(time_obs[i, o], log_y0_o, log_y1_o,
                                             t1_o, alpha_o, shape_o);
          log_lik[i] += normal_lpdf(log_y[i, o, k] | mu_log, tau_eps[k]);
        }
      }
    }
  }
}
