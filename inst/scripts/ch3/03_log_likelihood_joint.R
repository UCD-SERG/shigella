#!/usr/bin/env Rscript
# =====================================================================
# 03_log_likelihood_joint.R -- the shared-infection-time likelihood
#
# WHAT THIS IS FOR. serocalculator's log_likelihood() loops over biomarkers and
# adds the results, which per participant gives
#
#     [ int_t q1(y1|t) u(t) dt ] x [ int_t' q2(y2|t') u(t') dt' ]
#
# -- each measurement integrated over its own infection time. The methodology
# vignette states a different model (its equation for multiple biomarkers),
#
#     int_t q1(y1|t) q2(y2|t) u(t) dt
#
# and notes it "is not the product of the two single-biomarker likelihoods".
# This file implements that form, and pairs the Monte-Carlo draws across
# biomarkers, which relaxes the conditional-independence assumption of §9.2.
# That assumption holds exactly when the curve parameters are independent
# across biomarkers, which Chapter 2 shows they are not.
#
# THE NOISE MODEL, READ OFF src/serocalc.c RATHER THAN GUESSED
#
#   dnsB(y) = [prbF(y) - prbF(y-nu)]/nu
#       biological noise is ADDITIVE uniform:      z = y_f + B,  B ~ U(0, nu)
#   IdnsB(...)/(2 eps) integrates dnsB(z)/z over z
#       measurement noise is MULTIPLICATIVE uniform: y = z U, U ~ U(1-eps,1+eps)
#
#   Composing the two gives a closed form, so no inner quadrature is needed:
#
#       p(y | y_f) = log(zmax/zmin) / (2 eps nu)
#       zmin = max(y_f,      y/(1+eps))
#       zmax = min(y_f + nu, y/(1-eps))          zero when zmin >= zmax
#
#   Check: for a never-infected subject y_f = 0, which reduces to the
#   (Qa/nu) log(zmax/zmin)/(2 eps) term in dnsBM() exactly.
#
# LEFT CENSORING. serocalc.c uses prbBM() -- a probability, not a density --
# whenever y <= y.low. In this data 31% of HlyE IgG and 16% of HlyE IgA sit at
# the limit, so ignoring this is not a rounding matter. P(Y <= y.low | y_f) is
# computed by integrating the same closed form.
#
# UNITS. alpha is stored per day and age is in years, so alpha is multiplied by
# 365.25 before the quadrature, exactly as log_likelihood.R line 124 does.
#
# THE Pa / Qa MIXTURE. The density is a mixture over whether the participant has
# ever been infected, with weights Pa = 1 - exp(-lambda a) and Qa = exp(-lambda a).
# serocalc.c carries Pa inside dnsF and Qa inside the first branch of dnsBM.
# Omitting them leaves the mixture weights independent of lambda, which flattens
# the likelihood almost completely.
#
# THE `iter` COLUMN IS LOAD-BEARING. It is the only thing linking a person's
# IgG parameters to their IgA parameters. If it is missing or misaligned the
# pairing silently reverts to the product form and the chapter's contribution
# disappears with no error raised.
# =====================================================================

# ---------------------------------------------------------------------
# noise-free decay curve, serocalculator parameterisation
#   y_f(t) = y1 / (1 + d y1^d alpha t)^(1/d),   d = r - 1
# The rise is ignored, as the package does.
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# WHICH FORM.  Default reproduces the previous behaviour exactly.
#   product    A   separate tau, separate draws   (current practice)
#   sharedtau  B   one tau, draws averaged separately
#   joint      C   one tau, paired draws          (chapter 3)
# All three share every other component of this file, so a comparison
# between them cannot be contaminated by a convention mismatch.
# ---------------------------------------------------------------------
LIK_FORM <- Sys.getenv("LIK_FORM", "joint")
stopifnot(LIK_FORM %in% c("product", "sharedtau", "joint"))

.curve_mat <- function(t, y1, alpha, r) {
  d <- r - 1
  # outer product over draws x nodes, vectorised
  base <- 1 + outer(d * y1^d * alpha, t)
  y1 / base^(1 / d)
}

# ---------------------------------------------------------------------
# p(y | y_f), closed form. yf may be a matrix; y is scalar.
# ---------------------------------------------------------------------
.dens <- function(y, yf, nu, eps) {
  zmin <- pmax(yf, y / (1 + eps))
  zmax <- pmin(yf + nu, y / (1 - eps))
  ifelse(zmin < zmax, log(zmax / zmin) / (2 * eps * nu), 0)
}

# ---------------------------------------------------------------------
# P(Y <= y0 | y_f), for left-censored observations.
# Integrating the closed form has no elementary antiderivative in general, so
# this uses a short Gauss-Legendre rule over z; 24 nodes is ample because the
# integrand is smooth in z on a short interval.
# ---------------------------------------------------------------------
.GL <- local({
  n <- 24L
  # nodes and weights on (0,1)
  k <- seq_len(n - 1)
  b <- k / sqrt(4 * k^2 - 1)
  J <- diag(0, n)
  J[cbind(k, k + 1)] <- b
  J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = (rev(e$values) + 1) / 2, w = 2 * rev(e$vectors[1, ])^2 / 2)
})

.prob_below <- function(y0, yf, nu, eps) {
  # P(Y <= y0) = E_z[ P(U <= y0/z) ], z ~ U(y_f, y_f + nu), U ~ U(1-eps, 1+eps)
  out <- matrix(0, nrow(yf), ncol(yf))
  for (i in seq_along(.GL$x)) {
    z <- yf + nu * .GL$x[i]
    u <- pmin(1 + eps, pmax(1 - eps, y0 / z))
    out <- out + .GL$w[i] * (u - (1 - eps)) / (2 * eps)
  }
  out
}

# ---------------------------------------------------------------------
# quadrature grid on log(1 + t)
# Both the backward-recurrence density and the decay curve vary fastest near
# zero, so an even grid in t wastes most of its nodes.
# ---------------------------------------------------------------------
.grid <- function(a, n_tau = 400L) {
  if (n_tau %% 2 == 0) n_tau <- n_tau + 1L
  z <- seq(0, log1p(a), length.out = n_tau)
  h <- z[2] - z[1]
  w <- c(1, rep(c(4, 2), length.out = n_tau - 2), 1) * h / 3
  list(t = expm1(z), w = w * exp(z)) # exp(z) is dt/dz
}

# =====================================================================
# log_likelihood_joint()
# =====================================================================
log_likelihood_joint <- function(lambda, pop_data, curve_params, noise_params,
                                 antigen_isos = NULL, n_tau = 400L,
                                 verbose = FALSE) {
  # antigen_iso arrives as a factor from some sources and a character from
  # others; comparing two factors with different level sets raises an error
  # rather than returning FALSE, so everything is coerced on entry.
  pop_data <- as.data.frame(pop_data)
  curve_params <- as.data.frame(curve_params)
  noise_params <- as.data.frame(noise_params)
  pop_data$antigen_iso <- as.character(pop_data$antigen_iso)
  curve_params$antigen_iso <- as.character(curve_params$antigen_iso)
  noise_params$antigen_iso <- as.character(noise_params$antigen_iso)
  antigen_isos <- as.character(
    if (is.null(antigen_isos)) unique(pop_data$antigen_iso) else antigen_isos
  )
  K <- length(antigen_isos)

  if (!"iter" %in% names(curve_params)) {
    warning("curve_params has no `iter` column; the draws cannot be paired, ",
      "so this reduces to the product form.",
      call. = FALSE
    )
    curve_params$iter <- ave(seq_len(nrow(curve_params)),
      curve_params$antigen_iso,
      FUN = seq_along
    )
  }

  iters <- Reduce(
    intersect,
    lapply(antigen_isos, function(a) {
      curve_params$iter[curve_params$antigen_iso == a]
    })
  )
  # UNITS. The banks store alpha per DAY, which is how Chapter 2 estimates it,
  # while `age` is in YEARS and the quadrature grid runs from 0 to age. The
  # package handles this in log_likelihood.R line 124 with
  #     alpha = alpha * 365.25
  # and the same conversion has to happen here. Without it the curve decays 365
  # times too slowly, y_f never leaves the neighbourhood of the peak, and the
  # density is zero for almost every participant.
  par_list <- lapply(antigen_isos, function(a) {
    d <- curve_params[curve_params$antigen_iso == a, ]
    d <- d[match(iters, d$iter), c("y1", "alpha", "r")]
    d$alpha <- d$alpha * 365.25
    d
  })
  noise <- lapply(antigen_isos, function(a) {
    n <- noise_params[noise_params$antigen_iso == a, ]
    if (nrow(n) != 1) {
      stop("need exactly one noise row for ", a, " (got ",
        nrow(n), "); subset by Country first",
        call. = FALSE
      )
    }
    as.list(n)
  })
  M <- length(iters)
  if (verbose) message("  ", M, " paired draws, ", K, " biomarkers")

  # ---- the expensive part, done once ---------------------------------------
  # y_f depends only on the draw and the node, never on the person, so the
  # M x G matrices are built once and reused for everyone of the same age.
  # Ages are rounded to a grid so the cache is small; the curve is smooth in a.
  # Ages are binned before caching. The likelihood is smooth in age and the
  # curve matrix costs M x G doubles, so caching at 0.1 years would build a few
  # hundred matrices; 0.5 years cuts that fivefold at negligible cost.
  AGE_BIN <- as.numeric(Sys.getenv("AGE_BIN", "0.5"))
  # THE CENSORED BRANCH CACHES TOO, and that is where the time goes.
  # P(Y <= y.low) is a function of (y.low, y_f, nu, eps), and all four are fixed
  # once the age bin and the biomarker are fixed -- it does not depend on the
  # participant. Computing it per person meant running a 24-node quadrature over
  # a 1000 x 401 matrix for every censored observation, which is why Nepal (53%
  # censored) ran four times slower than Bangladesh (1.6%).
  cache <- new.env(parent = emptyenv())
  ccache <- new.env(parent = emptyenv())

  get_yf <- function(a) {
    key <- sprintf("%.2f", a)
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }
    g <- .grid(a, n_tau)
    yf <- lapply(seq_len(K), function(k) {
      .curve_mat(g$t, par_list[[k]]$y1, par_list[[k]]$alpha, par_list[[k]]$r)
    })
    val <- list(g = g, yf = yf, key = key)
    assign(key, val, envir = cache)
    val
  }

  get_cens <- function(cc, k) {
    key <- paste0(cc$key, "|", k)
    v <- ccache[[key]]
    if (!is.null(v)) {
      return(v)
    }
    nk <- noise[[k]]
    v <- .prob_below(nk$y.low, cc$yf[[k]], nk$nu, nk$eps)
    assign(key, v, envir = ccache)
    v
  }

  ids <- unique(pop_data$id)
  ll <- 0
  for (i in ids) {
    d <- pop_data[pop_data$id == i, ]
    y <- d$value[match(antigen_isos, d$antigen_iso)]
    if (anyNA(y)) next
    a <- round(d$age[1] / AGE_BIN) * AGE_BIN
    a <- max(a, AGE_BIN) # never zero: the grid needs a > 0
    cc <- get_yf(a)

    q <- vector("list", K)
    for (k in seq_len(K)) {
      nk <- noise[[k]]
      q[[k]] <- if (y[k] <= nk$y.low) {
        get_cens(cc, k) # left censored
      } else {
        .dens(y[k], cc$yf[[k]], nk$nu, nk$eps)
      }
    }
    # --- never-infected term: per biomarker, and as a product ----------
    never_k <- vapply(seq_len(K), function(k) {
      nk <- noise[[k]]
      if (y[k] <= nk$y.low) {
        .prob_below(nk$y.low, matrix(0, 1, 1), nk$nu, nk$eps)[1, 1]
      } else {
        .dens(y[k], 0, nk$nu, nk$eps)
      }
    }, numeric(1))
    never <- prod(never_k)

    # --- the Pa / Qa mixture, common to all three forms -----------------
    #   Pa = 1 - exp(-lambda a)   at least one infection by age a
    #   Qa = exp(-lambda a)       none; the noise-free value is then zero
    #   u  integrates to one over (0, a)
    Pa <- 1 - exp(-lambda * a)
    Qa <- exp(-lambda * a)
    u <- lambda * exp(-lambda * cc$g$t) + exp(-lambda * a) / a
    wu <- cc$g$w * u

    # --- and here, and ONLY here, the three forms differ ----------------
    f <- switch(LIK_FORM,

      # C: integrate the paired product for each draw, then average
      joint = mean(Pa * as.vector(Reduce(`*`, q) %*% wu)) + Qa * never,

      # B: average each biomarker over draws first, multiply, integrate.
      #    One tau is still shared; only the pairing is gone.
      sharedtau = Pa * sum(Reduce(`*`, lapply(q, colMeans)) * wu) +
        Qa * never,

      # A: each biomarker gets a complete single-marker density of its
      #    own -- its own tau, its own draws, its own mixture -- and the
      #    two are multiplied.  This is what summing log-likelihoods does.
      product = prod(vapply(
        seq_len(K), function(k) {
          Pa * sum(colMeans(q[[k]]) * wu) + Qa * never_k[k]
        },
        numeric(1)
      )),
      stop("LIK_FORM must be product, sharedtau or joint", call. = FALSE)
    )
    ll <- ll + log(max(f, .Machine$double.xmin))
  }
  ll
}

# =====================================================================
# THE GATE: with one biomarker, joint must equal the package's likelihood
# =====================================================================
validate_joint_likelihood <- function(pop_data, pop_pkg, curve_params,
                                      noise_params, antigen_isos,
                                      lambdas = c(0.02, 0.1, 0.5),
                                      n_tau = 400L, tol = 1e-3) {
  cat("\n  one biomarker: joint must equal serocalculator's log_likelihood\n\n")
  pop_data <- as.data.frame(pop_data)
  curve_params <- as.data.frame(curve_params)
  noise_params <- as.data.frame(noise_params)
  for (nm in c("pop_data", "curve_params", "noise_params")) {
    d <- get(nm)
    d$antigen_iso <- as.character(d$antigen_iso)
    assign(nm, d)
  }
  antigen_isos <- as.character(antigen_isos)

  res <- list()
  for (a in antigen_isos) {
    pd <- pop_data[pop_data$antigen_iso == a, ]
    pdk <- pop_pkg[as.character(pop_pkg$antigen_iso) == a, ]
    cp <- curve_params[curve_params$antigen_iso == a, ]
    np <- noise_params[noise_params$antigen_iso == a, ]
    for (lam in lambdas) {
      pkg <- serocalculator::log_likelihood(
        lam, pdk, serocalculator::as_sr_params(cp), np,
        antigen_isos = a
      )
      own <- log_likelihood_joint(lam, pd, cp, np, antigen_isos = a, n_tau = n_tau)
      d <- abs(pkg - own)
      pass <- is.finite(d) && d < tol * max(1, abs(pkg))
      cat(sprintf(
        "    %-10s lambda %5.3f   package %14.4f   joint %14.4f   diff %9.2e  %s\n",
        a, lam, pkg, own, d, if (pass) "ok" else "FAIL"
      ))
      res[[length(res) + 1]] <- data.frame(
        antigen_iso = a, lambda = lam,
        package = pkg, joint = own,
        diff = d, pass = pass
      )
    }
  }
  out <- do.call(rbind, res)
  cat(sprintf("\n    %d of %d checks pass\n\n", sum(out$pass), nrow(out)))
  out
}
