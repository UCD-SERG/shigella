#!/usr/bin/env Rscript
# =====================================================================
# 43_three_form_coverage.R  --  A, B and C side by side
#
#   cd ~/"chapter 2 work"/sees/ch3
#   Rscript R/43_three_form_coverage.R
#
# Reads out/cov (forms A and C) and out/covB (form B) and pairs them on
# (lambda, seed, replicate).  The three were run at the same seed, so
# replicate r is the same simulated dataset in all three and the
# comparison is paired rather than between independent samples.
#
# WHAT THE ANSWER MEANS
#   A is 0.897 and C is 0.980.  Where B lands says how the repair divides
#   between the two assumptions:
#
#     B close to C   sharing one infection time does the repair on its own,
#                    and pairing the draws adds width without adding
#                    calibration.  Chapter 2's contribution to chapter 3
#                    would then be smaller than the width contrast
#                    suggests, and the chapter should say so.
#
#     B between      the pairing supplies the remainder.  Chapter 2 earns
#                    its place on the evidence that matters.
#
#   Either is reportable.  Neither was known before this run.
# =====================================================================

suppressPackageStartupMessages(library(dplyr))
say <- function(...) cat(..., "\n", sep = "")
hr <- function(c = "-") cat(strrep(c, 76), "\n")
dir.create("tab", showWarnings = FALSE)

rd <- function(f) {
  r <- try(readRDS(f), silent = TRUE)
  if (inherits(r, "try-error")) NULL else as.data.frame(r)
}
load_dir <- function(d) {
  f <- list.files(d, "^cov_.*\\.rds$", full.names = TRUE)
  if (!length(f)) {
    return(NULL)
  }
  do.call(rbind, lapply(f, rd))
}

hr("=")
say("THREE FORMS, PAIRED ON THE SAME SIMULATED DATASETS")
hr("=")
AC <- load_dir("out/cov")
B <- load_dir("out/covB")
if (is.null(AC)) stop("no files in out/cov", call. = FALSE)
if (is.null(B)) say("  ⚠ out/covB is empty -- form B has not been run")

D <- rbind(AC, B)
D$form <- factor(
  ifelse(D$form == "composite", "A product",
    ifelse(D$form == "sharedtau", "B shared time", "C joint")
  ),
  levels = c("A product", "B shared time", "C joint")
)
D <- D[D$npop == 300, ]

# ---------------------------------------------------------------------
say(sprintf("  %d rows; replicates per (lambda, form):", nrow(D)))
print(table(lambda = D$lambda, form = D$form))

# only lambdas where all three are present can be compared
ok <- D %>%
  group_by(lambda) %>%
  summarise(k = n_distinct(form), .groups = "drop")
full <- ok$lambda[ok$k == 3]
say(
  "\n  lambdas with all three forms: ",
  if (length(full)) paste(full, collapse = ", ") else "none"
)

hr("=")
say("COVERAGE BY LAMBDA")
hr("=")
S <- D %>%
  group_by(lambda, form) %>%
  summarise(
    n = n(), coverage = mean(covers), width = mean(width),
    bias = 100 * (mean(est) / lambda[1] - 1), .groups = "drop"
  ) %>%
  mutate(mcse = sqrt(coverage * (1 - coverage) / n))
g <- function(l, f, v) {
  x <- S[[v]][S$lambda == l & S$form == f]
  if (length(x)) x[1] else NA_real_
}
lams <- sort(unique(S$lambda))
say(sprintf(
  "  %7s %10s %10s %10s   %10s %10s %10s", "lambda",
  "A", "B", "C", "w A", "w B", "w C"
))
for (l in lams) {
  say(sprintf(
    "  %7.2f %10.3f %10.3f %10.3f   %10.4f %10.4f %10.4f", l,
    g(l, "A product", "coverage"), g(l, "B shared time", "coverage"),
    g(l, "C joint", "coverage"),
    g(l, "A product", "width"), g(l, "B shared time", "width"),
    g(l, "C joint", "width")
  ))
}

# ---------------------------------------------------------------------
hr("=")
say("POOLED")
hr("=")
P <- D %>%
  group_by(form) %>%
  summarise(N = n(), coverage = mean(covers), width = mean(width), .groups = "drop") %>%
  mutate(
    mcse = sqrt(coverage * (1 - coverage) / N),
    z = (coverage - 0.95) / mcse
  )
say(sprintf(
  "  %-14s %7s %10s %8s %9s %10s", "form", "N", "coverage",
  "MCSE", "SE from .95", "mean width"
))
for (i in seq_len(nrow(P))) {
  with(
    P[i, ],
    say(sprintf("  %-14s %7d %10.3f %8.3f %+9.1f %10.4f", form, N, coverage, mcse, z, width))
  )
}
write.csv(P, "tab/table11_three_form_coverage.csv", row.names = FALSE)
write.csv(S, "tab/table11b_three_form_by_lambda.csv", row.names = FALSE)

# ---------------------------------------------------------------------
if ("B shared time" %in% P$form) {
  hr("=")
  say("HOW THE REPAIR DIVIDES")
  hr("=")
  a <- P$coverage[P$form == "A product"]
  b <- P$coverage[P$form == "B shared time"]
  cc <- P$coverage[P$form == "C joint"]
  say(sprintf("  A %.3f   ->   B %.3f   ->   C %.3f", a, b, cc))
  say(sprintf("\n  sharing one infection time   %+.3f", b - a))
  say(sprintf("  pairing the draws            %+.3f", cc - b))
  tot <- cc - a
  if (abs(tot) > 1e-9) {
    say(sprintf(
      "\n  sharing accounts for %.0f%% of the movement from A to C",
      100 * (b - a) / tot
    ))
  }
  wa <- P$width[P$form == "A product"]
  wb <- P$width[P$form == "B shared time"]
  wc <- P$width[P$form == "C joint"]
  say(sprintf(
    "\n  width   A %.4f  ->  B %.4f (%+.1f%%)  ->  C %.4f (%+.1f%%)",
    wa, wb, 100 * (wb / wa - 1), wc, 100 * (wc / wb - 1)
  ))
  cat("
  ★ Compare the coverage split with the width split reported on the
    observed data, which was about three quarters to sharing.  If the two
    agree, the attribution holds on simulated and observed data by two
    different measures.  If they do not, the widths and the calibration
    are telling different stories and the chapter should say which it
    trusts.

  ⚠ Read the width row as well as the coverage row.  A form can reach
    nominal coverage by being wide for the wrong reason, and B adding
    width without adding calibration would be exactly that.
")
} else {
  hr("=")
  say("form B not present -- run run_formB_coverage.sh first")
  hr("=")
}
hr("=")
say("  wrote tab/table11_three_form_coverage.csv and tab/table11b_three_form_by_lambda.csv")
