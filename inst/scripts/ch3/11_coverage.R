#!/usr/bin/env Rscript
# =====================================================================
# 11_coverage.R  --  the missing pillar of Chapter 3
#
# THE QUESTION
#   The joint interval is about 26% wider than the composite one.  Wider
#   is not the same as right.  The only way to show the joint width is
#   the CORRECT width is to build intervals on data where lambda is known
#   and count how often they contain it.
#
#   This matters more than usual because the conventional check failed:
#   the cluster-robust SE ratio averaged 0.99 over the six cells and fell
#   below one in three of them.  Coverage is the only remaining route.
#
# WHY A HAND-WRITTEN SIMULATOR
#   sim_pop_data() draws each antigen-isotype independently, which
#   destroys exactly the pairing under test.  The simulator here is the
#   estimator run backwards: draw tau, draw ONE bank index per person,
#   evaluate both curves at that tau, add the same two-stage noise the
#   likelihood assumes, censor at y.low.
#
# COST -- read this before launching anything
#   Measured: joint at n=30, N_LAM=11, N_TAU=200 costs 50 s on one core.
#   Everything scales linearly in n, in N_LAM and in N_TAU.  Composite is
#   about 1/35 of joint.  So:
#       one replicate, n=300, N_LAM=41, N_TAU=200   ~31 min
#       40 replicates in one process                ~21 h
#       40 replicates over 8 shards                 ~2.6 h wall
#   Shard it.  Do not run 40 replicates in one process.
#
# STAGES
#   selftest   composite only, n=2000.  ~12 min.  This is the real gate:
#              serocalculator's own likelihood must recover lambda from
#              our simulated data.  If it does, the simulator matches the
#              established model.
#   jointcheck joint only, n=150, one sample.  ~10 min.  Confirms the
#              joint implementation behaves on simulated data too.
#   run        replicates.  Use SHARD / NSHARD to parallelise.
#
#   STAGE=selftest   Rscript R/11_coverage.R
#   STAGE=jointcheck Rscript R/11_coverage.R
#   STAGE=run SHARD=1 NSHARD=8 REPS=40 Rscript R/11_coverage.R
# =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(serocalculator) })

DATA    <- Sys.getenv("DATA_DIR", "data")
OUT     <- Sys.getenv("OUT_DIR",  "out/cov")
BANK    <- Sys.getenv("BANK",     "arm_T")
COUNTRY <- Sys.getenv("COUNTRY",  "Bangladesh")
LAMBDA  <- as.numeric(Sys.getenv("LAMBDA", "0.1"))
NPOP    <- as.integer(Sys.getenv("NPOP",   "300"))
REPS    <- as.integer(Sys.getenv("REPS",   "40"))
N_TAU   <- as.integer(Sys.getenv("N_TAU",  "200"))
N_LAM   <- as.integer(Sys.getenv("N_LAM",  "41"))
SPAN    <- as.numeric(Sys.getenv("SPAN",   "2.5"))   # grid runs lambda/SPAN .. lambda*SPAN
SEED0   <- as.integer(Sys.getenv("SEED",   "1000"))
# TIME UNITS.  age and lambda are in YEARS; the bank's t1 and alpha are
# stored PER DAY (the joint likelihood multiplies alpha by 365.25 for
# exactly this reason -- that mismatch was error #3 of the four the K=1
# gate caught).  tau must therefore be converted before it reaches the
# curve.  Left as an environment variable so the self-test can settle it
# empirically instead of by argument.
TAU_SCALE <- as.numeric(Sys.getenv("TAU_SCALE", "365.25"))
STAGE   <- Sys.getenv("STAGE",   "selftest")
SHARD   <- as.integer(Sys.getenv("SHARD",  "1"))
NSHARD  <- as.integer(Sys.getenv("NSHARD", "1"))
ISOS    <- c("HlyE_IgG", "HlyE_IgA")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(..., "\n", sep = "")
hr  <- function(c = "-") cat(strrep(c, 74), "\n")

# The likelihood file reads LIK_FORM when it is sourced, so each form
# needs its own environment.  Built once and cached; sourcing on every
# call would dominate the run time.
.lik_env <- new.env(parent = emptyenv())
.lik <- function(form) {
  if (!is.null(.lik_env[[form]])) return(.lik_env[[form]])
  old <- Sys.getenv("LIK_FORM", unset = NA)
  Sys.setenv(LIK_FORM = form)
  e <- new.env(parent = globalenv())
  sys.source(file.path("R", "03_log_likelihood_joint.R"), envir = e)
  if (is.na(old)) Sys.unsetenv("LIK_FORM") else Sys.setenv(LIK_FORM = old)
  assign(form, e, envir = .lik_env)
  e
}
# sourced once at top level as well, so anything expecting the functions
# in the global environment still finds them
source(file.path("R", "03_log_likelihood_joint.R"))

# ---------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------
chr <- function(d) { d <- as.data.frame(d)
                     d$antigen_iso <- as.character(d$antigen_iso); d }

cp    <- chr(readRDS(file.path(DATA, sprintf("curve_params_%s.rds", BANK))))
xs    <- chr(readRDS(file.path(DATA, "sees_xs_subset.rds")))
noise <- chr(serocalculator::example_noise_params_sees) |>
           filter(Country == COUNTRY, antigen_iso %in% ISOS)
stopifnot(nrow(noise) == 2)

AGES <- unique(xs[, c("id", "age")])$age
AGES <- AGES[is.finite(AGES) & AGES > 0]

G <- cp[cp$antigen_iso == ISOS[1], ]; G <- G[order(G$iter), ]
A <- cp[cp$antigen_iso == ISOS[2], ]; A <- A[order(A$iter), ]
stopifnot(nrow(G) == nrow(A), all(G$iter == A$iter))
NDRAW <- nrow(G)

# grid: centred on the truth and only as wide as it needs to be.  A wide
# grid wastes resolution -- with SPAN 2.5 and 41 points the neighbours are
# 4.7% apart, matching the real analysis, over a range the interval
# cannot plausibly leave at these sample sizes.
LAMS  <- exp(seq(log(LAMBDA / SPAN), log(LAMBDA * SPAN), length.out = N_LAM))
STEPP <- 100 * (SPAN^(2 / (N_LAM - 1)) - 1)

hr("="); say("COVERAGE STUDY  --  stage: ", STAGE); hr("=")
say(sprintf("  bank %s   %d paired draws   columns: %s",
            BANK, NDRAW, paste(names(G), collapse = ", ")))
say(sprintf("  country %s   y.low %s   nu %s   eps %s", COUNTRY,
            paste(sprintf("%.3f", noise$y.low), collapse = "/"),
            paste(sprintf("%.3g", noise$nu),    collapse = "/"),
            paste(sprintf("%.3g", noise$eps),   collapse = "/")))
say(sprintf("  true lambda %.4f   grid %d points over %.4f to %.4f (%.1f%% apart)",
            LAMBDA, N_LAM, min(LAMS), max(LAMS), STEPP))
say(sprintf("  ages resampled from %d participants, median %.1f",
            length(AGES), median(AGES)))
say("  curve: DECAY-ONLY from the peak (y1, alpha, r), matching the likelihood")
say("  tau:   mixture  lambda*exp(-lambda*t) + exp(-lambda*a)/a  on (0, a)")
say(sprintf("  isotype states: %s",
            if (as.integer(Sys.getenv("TAU_INDEP", "0")) == 1L)
              "INDEPENDENT per isotype  (the product likelihood's own model)"
            else "shared within participant  (one infection, one tau)"))

# ---------------------------------------------------------------------
# DECAY-ONLY curve, matching .curve_mat() in 03_log_likelihood_joint.R and
# the three parameters serocalculator passes to C (y1, alpha, d = r - 1):
#
#     y_f(t) = y1 / (1 + d * y1^d * alpha * t)^(1/d)
#
# t is measured from the PEAK, not from infection.  The rise phase and the
# pre-infection level are not part of the estimator, so a simulator that
# includes them is generating from a different model -- which is exactly
# what the bias was.
#
# Computed on the log scale so y1^d cannot overflow.
# ---------------------------------------------------------------------
curve_ab <- function(t, y1, alpha, r) {
  d  <- r - 1
  lg <- log1p(exp(log(d) + log(alpha) + d * log(y1) + log(pmax(t, 0))))
  exp(log(y1) - lg / d)
}

# ---------------------------------------------------------------------
# one cross-sectional sample.  ONE bank index per person = the pairing.
# ---------------------------------------------------------------------
# TAU_INDEP=1 draws infection status, time since infection and the bank
# index SEPARATELY for each isotype, which is exactly the model the product
# likelihood assumes.  Under it the product form is correctly specified and
# must be unbiased; any bias that survives is a plumbing artefact rather
# than the cost of the shared-time assumption.  Default 0 = shared, the
# truth the rest of the chapter uses.
TAU_INDEP  <- as.integer(Sys.getenv("TAU_INDEP", "0"))

# the two suspect conventions.  Defaults = current behaviour.
NEVER_MODE <- Sys.getenv("NEVER_MODE", "zero")    # zero | y0
TAU_MODE   <- Sys.getenv("TAU_MODE",   "trexp")   # trexp | unif
stopifnot(NEVER_MODE %in% c("zero", "y0"), TAU_MODE %in% c("trexp", "unif"))

sim_xs <- function(lambda, n, seed) {
  set.seed(seed)
  age  <- sample(AGES, n, replace = TRUE)

  # one draw of (infected?, time since infection, bank index) for a person
  # u(tau) = lambda*exp(-lambda*tau) + exp(-lambda*a)/a   on (0, a)
  #
  # a mixture, not an atom plus a density: weight 1-exp(-lambda*a) on the
  # truncated exponential and weight exp(-lambda*a) on a uniform.  Everyone
  # sits on the curve; there is no never-infected branch.  This is the
  # convention serocalc.c uses and the one 03_log_likelihood_joint.R was
  # written against.
  # PA_QA -- the likelihood has TWO layers and both are reproduced here.
  #
  #   outer   Pa = 1 - exp(-lambda a)   infected, tau ~ u, on the curve
  #           Qa = exp(-lambda a)       never infected, y_f = 0
  #
  #   inner   u(t) = lambda*exp(-lambda t) + exp(-lambda a)/a   on (0, a)
  #           itself a mixture: truncated exponential with weight
  #           1-exp(-lambda a), uniform with weight exp(-lambda a)
  #
  # tau is returned in YEARS, because the likelihood builds its quadrature
  # grid on age in years and multiplies alpha by 365.25.  Matching that
  # explicitly is safer than using days against a per-day alpha and relying
  # on the product to come out equal.
  draw_state <- function() {
    Pa  <- 1 - exp(-lambda * age)
    ev  <- runif(n) < Pa                  # outer layer: infected at all?
    q   <- exp(-lambda * age)             # inner layer: weight on the uniform
    u   <- runif(n)
    uni <- runif(n) < q
    tt  <- ifelse(uni,
                  u * age,                                 # uniform on (0, a)
                  -log(1 - u * (1 - q)) / lambda)           # truncated exp
    list(ever = ev, tau_y = tt,           # YEARS
         k = sample.int(NDRAW, n, replace = TRUE))
  }

  st    <- draw_state()          # shared across isotypes unless TAU_INDEP
  ever  <- st$ever
  tau_y <- st$tau_y
  k     <- st$k

  out  <- vector("list", 2L)
  rise <- numeric(0)
  for (j in seq_along(ISOS)) {
    if (TAU_INDEP == 1L && j > 1L) {   # a fresh, unrelated state per isotype
      st <- draw_state(); ever <- st$ever; tau_y <- st$tau_y; k <- st$k
    }
    B  <- if (j == 1L) G else A
    nz <- noise[noise$antigen_iso == ISOS[j], ]
    # alpha is stored per DAY; the likelihood converts it to per YEAR and
    # evaluates the curve at tau in years.  Same conversion here.
    y <- rep(0, n)                                   # Qa branch: y_f = 0
    if (any(ever))
      y[ever] <- curve_ab(tau_y[ever], B$y1[k[ever]],
                          B$alpha[k[ever]] * 365.25, B$r[k[ever]])
    rise <- c(rise, mean(ever))          # now: share on the curve, not at zero
    y <- (y + runif(n, 0, nz$nu)) * runif(n, 1 - nz$eps, 1 + nz$eps)
    y <- pmax(y, nz$y.low)
    out[[j]] <- data.frame(id = seq_len(n), age = age, antigen_iso = ISOS[j],
                           value = y, Country = COUNTRY, llod = nz$y.low,
                           stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, out)
  attr(d, "frac_ever")  <- mean(ever)
  attr(d, "frac_rise")  <- mean(rise)
  attr(d, "med_tau")    <- if (any(ever)) median(tau_y[ever]) else NA_real_
  d
}

profile_ci <- function(f) {
  v <- vapply(LAMS, f, numeric(1)); i <- which.max(v)
  lo <- suppressWarnings(max(LAMS[seq_len(i)][v[seq_len(i)] < max(v) - 1.92]))
  hi <- suppressWarnings(min(LAMS[i:length(LAMS)][v[i:length(LAMS)] < max(v) - 1.92]))
  c(est = LAMS[i], lo = lo, hi = hi,
    edge = as.numeric(i == 1L || i == length(LAMS)))
}

# --- the three forms, all from 03_log_likelihood_joint.R --------------
.fit_form <- function(d, form) {
  e <- .lik(form)
  profile_ci(function(l) e$log_likelihood_joint(
    l, d, cp, noise, antigen_isos = ISOS, n_tau = N_TAU))
}
fit_product   <- function(d) .fit_form(d, "product")     # A
fit_sharedtau <- function(d) .fit_form(d, "sharedtau")   # B
fit_joint     <- function(d) .fit_form(d, "joint")       # C

# --- serocalculator, kept as a REFERENCE ONLY -------------------------
# Useful for showing that our machinery reproduces the established
# calculation at K = 1.  It is deliberately not the comparator: it differs
# from our code in the noise convention, the tau density and the
# never-infected branch, so a difference against it is not attributable
# to the independence assumption, which is the only thing being measured.
fit_pkg <- function(d) {
  pd  <- as_pop_data(d, age = "age", value = "value", id = "id",
                     antigen_isos = ISOS)
  srp <- as_sr_params(cp |> filter(antigen_iso %in% ISOS))
  profile_ci(function(l) serocalculator::log_likelihood(
    l, pd, srp, noise, antigen_isos = ISOS))
}

# `composite` now means OUR product form, not the package.  Anything that
# still asks for fit_composite() gets form A.
fit_composite <- fit_product

report <- function(nm, z) say(sprintf(
  "  %-10s est %.4f  [%.4f, %.4f]   truth %.4f   %s",
  nm, z["est"], z["lo"], z["hi"], LAMBDA,
  if (z["lo"] <= LAMBDA && LAMBDA <= z["hi"]) "covers" else "*** MISSES ***"))

# =====================================================================
# discord -- is the real data more discordant than arm_T predicts?
#
# THE OPEN ITEM.  Simulation says the composite OVERestimates lambda by
# about 17%.  The real data show the joint estimate ABOVE the composite by
# about 25%.  The signs disagree, and until that is explained no statement
# about the level of incidence can be made.
#
# One hypothesis explains it, together with two other loose ends -- the
# constant max_ll offset between banks, and the width ratio being larger in
# the real data (1.42) than in simulation (1.09-1.21):
#
#     real IgG-IgA pairs are MORE DISCORDANT than the chapter 2
#     correlation predicts.
#
# If true, forcing one tau onto a discordant pair pushes that tau toward
# the recent end, because early after infection both markers are high and
# can accommodate a wide range of value combinations, whereas an old tau
# requires both to be low and to have decayed together.  Simulated data are
# generated FROM arm_T, so they are exactly as concordant as arm_T says and
# the push is absent.
#
# The test needs no per-person tau: compare the observed joint distribution
# of the two isotypes against simulated data generated at that cell own
# fitted lambda, with the same n, the same ages and the same noise.
# =====================================================================
if (STAGE == "discord") {

  hr("="); say("DISCORDANCE -- observed versus simulated under arm_T"); hr("=")

  CELLS <- Sys.getenv("CELLS", paste(
    "highE|Bangladesh|0.3786", "highE|Nepal|0.0497", "highE|Pakistan|0.0922",
    "lowE|Bangladesh|0.2132",  "lowE|Nepal|0.0435",  "lowE|Pakistan|0.0567",
    sep = ","))
  NREP <- as.integer(Sys.getenv("DISC_REPS", "20"))

  # four statistics, none of which needs tau
  disc_stats <- function(d, lo) {
    g <- d$value[d$antigen_iso == ISOS[1]]
    a <- d$value[d$antigen_iso == ISOS[2]]
    cg <- g <= lo[1] + 1e-12; ca <- a <= lo[2] + 1e-12
    ok <- !cg & !ca
    r  <- if (sum(ok) > 10) cor(log(g[ok]), log(a[ok])) else NA_real_
    rs <- if (sum(ok) > 10) {
      f <- lm(log(a[ok]) ~ log(g[ok])); sd(residuals(f))
    } else NA_real_
    # "crossed": one isotype above its own median, the other below
    hi_g <- g > median(g); hi_a <- a > median(a)
    c(cor = r, resid_sd = rs, both_at_limit = mean(cg & ca),
      crossed = mean(xor(hi_g, hi_a)))
  }

  say(sprintf("\n  %-22s %6s %8s %9s %9s %9s %9s",
              "cell", "n", "source", "cor(gA)", "resid SD", "both low", "crossed"))

  for (cell in strsplit(CELLS, ",")[[1]]) {
    p <- strsplit(cell, "\\|")[[1]]
    arm <- p[1]; ctry <- p[2]; lam <- as.numeric(p[3])

    real <- xs[xs$arm == arm & xs$Country == ctry, ]
    if (!nrow(real)) { say("  skip ", cell, ": no rows"); next }

    # the simulator reads these globals; point them at this cell
    noise <<- chr(serocalculator::example_noise_params_sees) |>
                dplyr::filter(Country == ctry, antigen_iso %in% ISOS)
    lo <- c(noise$y.low[noise$antigen_iso == ISOS[1]],
            noise$y.low[noise$antigen_iso == ISOS[2]])
    AGES <<- unique(real[, c("id", "age")])$age
    n_i  <- length(unique(real$id))

    o <- disc_stats(real, lo)
    say(sprintf("  %-22s %6d %8s %9.3f %9.3f %8.1f%% %8.1f%%",
                paste(arm, ctry), n_i, "observed",
                o["cor"], o["resid_sd"], 100*o["both_at_limit"], 100*o["crossed"]))

    S <- vapply(seq_len(NREP), function(b)
      disc_stats(sim_xs(lam, n_i, SEED0 + 7000 + b), lo), numeric(4))
    m <- rowMeans(S, na.rm = TRUE); s <- apply(S, 1, sd, na.rm = TRUE) / sqrt(NREP)
    say(sprintf("  %-22s %6s %8s %9.3f %9.3f %8.1f%% %8.1f%%   (MCSE %.3f)",
                "", "", "arm_T sim", m[1], m[2], 100*m[3], 100*m[4], s[1]))
    say(sprintf("  %-22s %6s %8s %+9.3f %+9.3f %+8.1f%% %+8.1f%%   <- observed minus simulated",
                "", "", "difference",
                o["cor"]-m[1], o["resid_sd"]-m[2],
                100*(o["both_at_limit"]-m[3]), 100*(o["crossed"]-m[4])))
    cat("\n")
  }

  cat("
  HOW TO READ IT

    observed cor LOWER and resid SD HIGHER than simulated, and crossed HIGHER
        -> the hypothesis holds.  Real pairs are more discordant than arm_T
           predicts, which explains the sign disagreement, the constant
           max_ll offset and the larger real-data width ratio, all at once.

    the two agree
        -> a different explanation is needed.  Candidates, in order:
           repeated infection raising the baseline; age-varying kinetics;
           the bank coming from a case cohort rather than the general
           population.

    Note that cor() here is inflated by lambda itself -- at higher incidence
    more people are infected and both markers are elevated together.  That
    is why each cell is simulated at ITS OWN fitted lambda rather than at a
    common value, so the comparison is like for like.
")
  quit(status = 0)
}

# =====================================================================
# bigN -- finite-sample bias, or inconsistency?
#
# The composite bias barely moved when n doubled: +16.9% at n = 300 against
# +15.5% at n = 600, a ratio of 0.92 where 1/n would give 0.50 and 1/sqrt(n)
# 0.71.  That is the signature of an estimator converging to the wrong
# value rather than of a bias that vanishes with more data.
#
# The distinction changes the claim.  "Biased at this sample size" is a
# caveat; "does not go away however large the study" is a finding.
#
# Composite-only, so a very large n is cheap -- the composite costs about a
# thirty-fifth of the joint, and cost is linear in n.
# =====================================================================
if (STAGE == "bigN") {

  # which likelihood.  Default reproduces the previous behaviour exactly.
  FORM1 <- Sys.getenv("FORM1", "composite")
  if (!FORM1 %in% c("composite", "joint"))
    stop("FORM1 must be composite or joint", call. = FALSE)

  hr("="); say("LARGE-SAMPLE BIAS -- ", toupper(FORM1), " LIKELIHOOD"); hr("=")
  say(sprintf("  lambda %.4f   n %d   %d replicates   %s only",
              LAMBDA, NPOP, REPS, FORM1))
  say(sprintf("  grid %d points over %.4f to %.4f (%.2f%% apart)",
              N_LAM, min(LAMS), max(LAMS), STEPP))
  cat("\n")

  est <- numeric(REPS); wid <- numeric(REPS); cov <- integer(REPS)
  for (b in seq_len(REPS)) {
    t0 <- Sys.time()
    d  <- sim_xs(LAMBDA, NPOP, SEED0 + 9000 + b)
    z  <- if (FORM1 == "joint") fit_joint(d) else fit_composite(d)
    est[b] <- z["est"]; wid[b] <- z["hi"] - z["lo"]
    cov[b] <- as.integer(z["lo"] <= LAMBDA && LAMBDA <= z["hi"])
    say(sprintf("  rep %2d  est %.4f  bias %+6.1f%%  width %.4f  %s  %.1f min",
                b, z["est"], 100*(z["est"]/LAMBDA - 1), wid[b],
                if (cov[b] == 1) "covers" else "MISSES",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  bias <- mean(est)/LAMBDA - 1
  se   <- sd(est)/sqrt(REPS)/LAMBDA
  hr("="); say("RESULT"); hr("=")
  say(sprintf("  n %d   mean estimate %.4f   bias %+.1f%%  (MCSE %.1f%%)  = %.1f SE",
              NPOP, mean(est), 100*bias, 100*se, abs(bias)/se))
  say(sprintf("  mean width %.4f   coverage %d/%d", mean(wid), sum(cov), REPS))
  cat("
  COMPARE AGAINST THE SMALL-SAMPLE RUNS

    n =  300   bias +16.9%
    n =  600   bias +15.5%

    bias holds near +15%   -> INCONSISTENT.  The composite converges to the
                              wrong value and no sample size fixes it.  The
                              chapter can say so.
    bias falls toward 0    -> finite-sample bias.  The claim becomes
                              \"biased at realistic sample sizes\", which is
                              still worth reporting but is a smaller claim.

    A bias that is many MCSE from zero at this n settles it either way.
")
  saveRDS(data.frame(lambda = LAMBDA, npop = NPOP, rep = seq_len(REPS),
                     est = est, width = wid, covers = cov),
          file.path(OUT, sprintf("bigN_%s_lam%s_n%d.rds",
                                 FORM1, sub("\\.", "p", LAMBDA), NPOP)))
  quit(status = 0)
}

# =====================================================================
# selftest -- composite only, large n.  THE GATE.
# =====================================================================
if (STAGE == "selftest") {
  hr("="); say("SELF TEST -- form C must recover a lambda drawn from its own model"); hr("=")
  say("  The comparison in this chapter is internal: forms A, B and C come")
  say("  from one function and differ only in the order of sum, integral and")
  say("  product.  So the question is no longer whether our simulator")
  say("  reproduces serocalculator -- it is whether the simulator and form C")
  say("  are the same model stated in two directions.")
  say("  serocalculator is evaluated alongside as a reference, and a gap")
  say("  against it is expected rather than disqualifying.")
  cat("
  THE THING BEING SETTLED HERE
    age and lambda are in YEARS.  The bank's t1 and alpha are stored PER
    DAY -- the joint likelihood multiplies alpha by 365.25 for exactly
    this reason, and getting it wrong was error #3 of the four the K=1
    gate caught.  Rather than argue about which convention the curve
    wants, run both and see which one recovers the truth.

    TAU_SCALE = 365.25  tau in days   (expected correct)
    TAU_SCALE = 1       tau in years

    (The old rise-branch tell-tale no longer applies: the estimator has
     no rise phase, so the simulator has none either.  What to watch now
     is simply whether the ratio sits inside 3% at all three lambdas.)

")
  scales <- as.numeric(strsplit(Sys.getenv("SCALES", "365.25,1"), ",")[[1]])
  lams   <- as.numeric(strsplit(Sys.getenv("LAMS_TEST", "0.05,0.1,0.2"), ",")[[1]])
  N_SELF <- as.integer(Sys.getenv("N_SELF", "2000"))

  say(sprintf("  %9s %8s %9s %8s %9s %11s %8s", "TAU_SCALE", "lambda",
              "est", "ratio", "on rise", "at limit G", "min"))
  best <- NULL
  for (sc in scales) {
    TAU_SCALE <<- sc
    ratios <- numeric(0)
    for (lam in lams) {
      LAMS  <<- exp(seq(log(lam / SPAN), log(lam * SPAN), length.out = N_LAM))
      t0 <- Sys.time()
      d  <- sim_xs(lam, N_SELF, SEED0)
      z  <- fit_joint(d)                       # form C -- the gate
      zp <- try(fit_pkg(d), silent = TRUE)     # package -- reference only
      if (!inherits(zp, "try-error"))
        say(sprintf("        (reference: serocalculator %.4f, ratio %.2f)",
                    zp["est"], zp["est"] / lam))
      cg <- mean(d$value[d$antigen_iso == ISOS[1]] <=
                 noise$y.low[noise$antigen_iso == ISOS[1]] + 1e-12)
      ratios <- c(ratios, z["est"] / lam)
      say(sprintf("  %9.2f %8.3f %9.4f %8.2f %8.0f%% %10.1f%% %8.1f",
                  sc, lam, z["est"], z["est"] / lam,
                  100 * attr(d, "frac_rise"), 100 * cg,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    dev <- max(abs(log(ratios)))
    if (is.null(best) || dev < best$dev) best <- list(sc = sc, dev = dev,
                                                      ratios = ratios)
    cat("\n")
  }

  hr("="); say("VERDICT"); hr("=")
  say(sprintf("  best convention: TAU_SCALE = %g   ratios %s",
              best$sc, paste(sprintf("%.2f", best$ratios), collapse = ", ")))

  # A ratio equal to SPAN, or to 1/SPAN, is the grid boundary rather than an
  # estimate.  This morning every cell returned exactly lambda*SPAN and the
  # verdict still read like an estimate.
  onedge <- any(abs(best$ratios - SPAN) < 1e-6) ||
            any(abs(best$ratios - 1 / SPAN) < 1e-6)
  if (onedge) {
    say(sprintf("\n  *** ON THE GRID BOUNDARY (SPAN = %g).  Not an estimate. ***", SPAN))
    say("  The likelihood is pushing lambda past the edge of the grid, which")
    say("  means the simulated data do not look like anything the estimator")
    say("  expects.  Widen SPAN only to confirm it is not a grid artefact;")
    say("  the cause is a model mismatch, not the grid.")
    unlink(".selftest_pass")
    quit(status = 1)
  }

  unlink(".selftest_pass")
  if (best$dev < log(1.03)) {
    say("  PASS -- the estimate tracks the truth within 3% at every lambda.")
    writeLines(format(Sys.time()), ".selftest_pass")
    say(sprintf("  Put TAU_SCALE=%g in the run and grid stages.", best$sc))
  } else {
    say("  FAIL -- neither convention recovers the truth.")
    cat("
  Remaining suspects, in order, all of them local and all with the answer
  written in R/03_log_likelihood_joint.R:

    1. noise composition.  This simulator does (curve + U(0, nu)) * U(1-eps, 1+eps).
       The closed form log(zmax/zmin)/(2*eps*nu) implies exactly that order;
       if the likelihood does it the other way round, swap the two lines.
    2. censoring.  This simulator does pmax(y, y.low).  If the likelihood
       instead drops or reassigns censored values, match it.
    3. the tau mixture.  This simulator uses P(never) = exp(-lambda*age)
       and a truncated exponential otherwise.  Check the likelihood's
       age-conditional density for tau.

  A ratio that is CONSTANT across the three lambdas points at a scale
  error; a ratio that DRIFTS with lambda points at the tau mixture.
")
  }
  quit(status = 0)
}

# =====================================================================
# jointcheck -- joint only, small n
# =====================================================================
if (STAGE == "jointcheck") {
  hr("="); say("JOINT CHECK -- joint likelihood, n = 150, one sample"); hr("=")
  d <- sim_xs(LAMBDA, 150L, SEED0)
  t0 <- Sys.time()
  report("composite", fit_composite(d))
  report("joint",     fit_joint(d))
  say(sprintf("\n  %.1f min for n=150; a replicate at n=%d will cost about %.0f min",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), NPOP,
              as.numeric(difftime(Sys.time(), t0, units = "mins")) * NPOP / 150))
  cat("
  At n = 150 a single interval is wide and one sample proves nothing
  about coverage.  What this checks is that both estimators run on
  simulated data and give sane numbers, and it measures the per-replicate
  cost so the shard count can be chosen properly.
")
  quit(status = 0)
}

# =====================================================================
# run -- replicates, shardable
# =====================================================================
myreps <- seq_len(REPS)[(seq_len(REPS) - 1) %% NSHARD + 1 == SHARD]
hr("="); say(sprintf("REPLICATES  shard %d of %d  ->  %d reps: %s",
                     SHARD, NSHARD, length(myreps),
                     paste(range(myreps), collapse = "-"))); hr("=")

# the shard file is rewritten after EVERY replicate.  A shard holds
# 12-13 replicates at about 13 min each, so finishing before the first
# save would mean risking three hours of work on one process staying
# alive.  It also lets the report be run on partial results.
# The form set goes in the name only when it is not the default pair, so
# every file written so far keeps the name the report already expects.
#
# Braces are not optional here.  At top level R completes `if (cond) ""`
# at the end of the line and the `else` on the next one has nothing to
# attach to -- which is exactly how the first version of this patch
# failed.
.ftag <- Sys.getenv("FORMS_RUN", "composite,joint")
.ftag <- if (identical(.ftag, "composite,joint")) {
  ""
} else {
  paste0("_", gsub(",", "-", .ftag))
}
OUTF <- file.path(OUT, sprintf("cov_%s%s_lam%s_n%d_sd%d_s%02d.rds",
                               BANK, .ftag, sub("\\.", "p", LAMBDA),
                               NPOP, SEED0, SHARD))

rows <- list()
for (rep in myreps) {
  t0 <- Sys.time()
  d  <- sim_xs(LAMBDA, NPOP, SEED0 + rep)
  # FORMS_RUN picks which of the three to fit.  Default is the pair the
  # chapter compares; "sharedtau" alone reproduces the same datasets at
  # the same seed, so it can be merged with an existing A/C run.
  .want <- strsplit(Sys.getenv("FORMS_RUN", "composite,joint"), ",")[[1]]
  .all  <- list(composite = fit_product, sharedtau = fit_sharedtau,
                joint     = fit_joint)
  if (!all(.want %in% names(.all)))
    stop("FORMS_RUN must be composite, sharedtau or joint", call. = FALSE)
  f <- lapply(.all[.want], function(fn) fn(d))
  for (nm in names(f)) {
    z <- f[[nm]]
    rows[[length(rows) + 1]] <- data.frame(
      rep = rep, form = nm, est = z["est"], lo = z["lo"], hi = z["hi"],
      width = z["hi"] - z["lo"], edge = z["edge"],
      covers = as.integer(z["lo"] <= LAMBDA && LAMBDA <= z["hi"]),
      row.names = NULL)
  }
  say(sprintf("  rep %3d  %s  %.0f min", rep,
              paste(vapply(names(f), function(nm) sprintf(
                "%s %.4f [%.4f,%.4f] %s", substr(nm, 1, 4), f[[nm]]["est"],
                f[[nm]]["lo"], f[[nm]]["hi"],
                if (f[[nm]]["lo"] <= LAMBDA && LAMBDA <= f[[nm]]["hi"]) "Y" else "n"),
                character(1)), collapse = "   "),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  part <- do.call(rbind, rows)
  part$lambda <- LAMBDA; part$npop <- NPOP; part$seed0 <- SEED0
  saveRDS(part, OUTF)
}

res <- do.call(rbind, rows)
res$lambda <- LAMBDA; res$npop <- NPOP; res$seed0 <- SEED0
saveRDS(res, OUTF)
say("\n  saved -> ", OUTF, "  (", nrow(res) / 2, " replicates)")
say("  combine every shard with:  Rscript R/12_coverage_report.R")
