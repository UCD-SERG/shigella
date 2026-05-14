// ============================================================
// model_1.stan
//
// Model 1 (formerly Model A): Independent biomarker antibody kinetics.
// Each biomarker fit independently — no cross-biomarker correlation.
//
// Compatible with serodynamics::prep_data_stan() output:
//   N, K, P, max_obs, n_obs[N], time_obs[N, max_obs], log_y[N, max_obs, K]
//
// Compatible with prep_priors_stan(model = "model_1") output.
// ============================================================

functions {
  real two_phase_curve(real t, real y0, real y1, real t1,
                       real alpha, real shape) {
    real beta;
    if (t <= t1) {
      beta = log(y1 / y0) / t1;
      return y0 * exp(beta * t);
    } else {
      real base = pow(y1, 1 - shape) - (1 - shape) * alpha * (t - t1);
      if (base <= 0) {
        return y0 * 0.01;
      }
      return pow(base, 1 / (1 - shape));
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
  real<lower=0> lkj_eps_eta;
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
          real y0_o = exp(theta[i, k][1]);
          real y1_o = y0_o + exp(theta[i, k][2]);
          real t1_o = exp(theta[i, k][3]);
          real alpha_o = exp(theta[i, k][4]);
          real shape_o = exp(theta[i, k][5]) + 1;
          real mu_log = log(two_phase_curve(time_obs[i, o], y0_o, y1_o,
                                             t1_o, alpha_o, shape_o));
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
          real y0_o = exp(theta[i, k][1]);
          real y1_o = y0_o + exp(theta[i, k][2]);
          real t1_o = exp(theta[i, k][3]);
          real alpha_o = exp(theta[i, k][4]);
          real shape_o = exp(theta[i, k][5]) + 1;
          real mu_log = log(two_phase_curve(time_obs[i, o], y0_o, y1_o,
                                             t1_o, alpha_o, shape_o));
          log_lik[i] += normal_lpdf(log_y[i, o, k] | mu_log, tau_eps[k]);
        }
      }
    }
  }
}
