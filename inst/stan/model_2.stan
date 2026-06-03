// ============================================================
// model_2.stan — JAGS-ALIGNED version (v3)
//
// Key change from previous version: compute log(y(t)) DIRECTLY
// (matching the JAGS reference model.jags from Chapter 1),
// rather than computing y(t) then taking log().
//
// This avoids the exp -> arithmetic -> log round-trip, which:
//   1. Avoids overflow when y1 is large
//   2. Removes the need for a discontinuous fallback (no `if base <= 0`)
//   3. Aligns numerically with the JAGS model that works for Chapter 1
//
// Also: y1 = y0 + exp(Theta[k,2]), so log(y1) = log_sum_exp(log_y0, Theta[k,2])
// This is more stable than log(exp(.) + exp(.)).
//
// Model 2: Kronecker correlated antibody kinetics.
// vec(Theta_i) ~ MVN_KP(vec(M), Sigma_B kron Sigma_P)
// log y[i,o,1:K] ~ MVN_K(mu_log[i,o,1:K], Sigma_eps)
// ============================================================

functions {
  // Compute log(y(t)) DIRECTLY (matching JAGS reference)
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

  matrix kron_chol(matrix L_B, matrix L_P) {
    int K = rows(L_B);
    int P = rows(L_P);
    matrix[K * P, K * P] L_out = rep_matrix(0, K * P, K * P);
    for (i in 1:K) {
      for (j in 1:i) {
        for (p in 1:P) {
          for (q in 1:p) {
            L_out[(i - 1) * P + p, (j - 1) * P + q] = L_B[i, j] * L_P[p, q];
          }
        }
      }
    }
    return L_out;
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
  real<lower=0> tau_B_scale;
  real<lower=0> tau_eps_scale;
  real<lower=0> lkj_P_eta;
  real<lower=0> lkj_B_eta;
  real<lower=0> lkj_eps_eta;
}

parameters {
  matrix[K, P] M;

  cholesky_factor_corr[K] L_Omega_B;
  cholesky_factor_corr[P] L_Omega_P;
  vector<lower=0>[K] tau_B;
  vector<lower=0>[P] tau_P;

  cholesky_factor_corr[K] L_Omega_eps;
  vector<lower=0>[K] tau_eps;

  matrix[N, K * P] Z;
}

transformed parameters {
  array[N] matrix[K, P] Theta;

  matrix[K, K] L_Sigma_B = diag_pre_multiply(tau_B, L_Omega_B);
  matrix[P, P] L_Sigma_P = diag_pre_multiply(tau_P, L_Omega_P);
  matrix[K * P, K * P] L_kron = kron_chol(L_Sigma_B, L_Sigma_P);

  vector[K * P] mu_vec;
  for (k in 1:K) {
    for (p in 1:P) {
      mu_vec[(k - 1) * P + p] = M[k, p];
    }
  }

  for (i in 1:N) {
    vector[K * P] theta_vec = mu_vec + L_kron * to_vector(Z[i]);
    for (k in 1:K) {
      for (p in 1:P) {
        Theta[i, k, p] = theta_vec[(k - 1) * P + p];
      }
    }
  }
}

model {
  // Priors on M (matches JAGS mu.par ~ dmnorm structure)
  for (k in 1:K) {
    for (p in 1:P) {
      M[k, p] ~ normal(mu_hyp_mean[p], mu_hyp_sd[p]);
    }
  }

  L_Omega_B ~ lkj_corr_cholesky(lkj_B_eta);
  L_Omega_P ~ lkj_corr_cholesky(lkj_P_eta);
  L_Omega_eps ~ lkj_corr_cholesky(lkj_eps_eta);

  tau_B ~ cauchy(0, tau_B_scale);
  tau_P ~ cauchy(0, tau_P_scale);
  tau_eps ~ cauchy(0, tau_eps_scale);

  to_vector(Z) ~ std_normal();

  matrix[K, K] L_Sigma_eps = diag_pre_multiply(tau_eps, L_Omega_eps);

  // Likelihood — compute log(y(t)) directly (matching JAGS)
  for (i in 1:N) {
    if (n_obs[i] > 0) {
      for (o in 1:n_obs[i]) {
        vector[K] mu_log_o;
        vector[K] y_log_o;
        for (k in 1:K) {
          // Theta[i, k, *] holds:
          //   [1] = log(y0)
          //   [2] = log(y1 - y0)        -> y1 = y0 + exp(par[2])
          //   [3] = log(t1)
          //   [4] = log(alpha)
          //   [5] = log(shape - 1)      -> shape = exp(par[5]) + 1
          real log_y0_o = Theta[i, k, 1];
          // y1 = y0 + exp(par2) = exp(log_y0) + exp(par2)
          // log(y1) = log_sum_exp(log_y0, par2)  -- numerically stable
          real log_y1_o = log_sum_exp(log_y0_o, Theta[i, k, 2]);
          real t1_o     = exp(Theta[i, k, 3]);
          real alpha_o  = exp(Theta[i, k, 4]);
          real shape_o  = exp(Theta[i, k, 5]) + 1;

          mu_log_o[k] = log_two_phase_curve(time_obs[i, o],
                                             log_y0_o, log_y1_o,
                                             t1_o, alpha_o, shape_o);
          y_log_o[k] = log_y[i, o, k];
        }
        y_log_o ~ multi_normal_cholesky(mu_log_o, L_Sigma_eps);
      }
    }
  }
}

generated quantities {
  // NOTE: Kinetics recomputed here intentionally — Stan scoping requires
  // local variables to be redefined; this is not duplication that can be eliminated.
  corr_matrix[K] Omega_B = multiply_lower_tri_self_transpose(L_Omega_B);
  corr_matrix[P] Omega_P = multiply_lower_tri_self_transpose(L_Omega_P);
  corr_matrix[K] Omega_eps = multiply_lower_tri_self_transpose(L_Omega_eps);
  cov_matrix[K] Sigma_B = quad_form_diag(Omega_B, tau_B);
  cov_matrix[P] Sigma_P = quad_form_diag(Omega_P, tau_P);
  cov_matrix[K] Sigma_eps = quad_form_diag(Omega_eps, tau_eps);

  array[N, K] real y0;
  array[N, K] real y1;
  array[N, K] real t1;
  array[N, K] real alpha;
  array[N, K] real shape;
  for (i in 1:N) {
    for (k in 1:K) {
      y0[i, k]    = exp(Theta[i, k, 1]);
      y1[i, k]    = y0[i, k] + exp(Theta[i, k, 2]);
      t1[i, k]    = exp(Theta[i, k, 3]);
      alpha[i, k] = exp(Theta[i, k, 4]);
      shape[i, k] = exp(Theta[i, k, 5]) + 1;
    }
  }

  vector[N] log_lik;
  {
    // L_Sigma_eps rebuilt here per Stan scoping; intentional, not accidental.
    matrix[K, K] L_Sigma_eps = diag_pre_multiply(tau_eps, L_Omega_eps);
    for (i in 1:N) {
      log_lik[i] = 0;
      if (n_obs[i] > 0) {
        for (o in 1:n_obs[i]) {
          vector[K] mu_log_o;
          vector[K] y_log_o;
          for (k in 1:K) {
            real log_y0_o = Theta[i, k, 1];
            real log_y1_o = log_sum_exp(log_y0_o, Theta[i, k, 2]);
            real t1_o     = exp(Theta[i, k, 3]);
            real alpha_o  = exp(Theta[i, k, 4]);
            real shape_o  = exp(Theta[i, k, 5]) + 1;
            mu_log_o[k] = log_two_phase_curve(time_obs[i, o],
                                               log_y0_o, log_y1_o,
                                               t1_o, alpha_o, shape_o);
            y_log_o[k] = log_y[i, o, k];
          }
          log_lik[i] += multi_normal_cholesky_lpdf(y_log_o | mu_log_o, L_Sigma_eps);
        }
      }
    }
  }
}
