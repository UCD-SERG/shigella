#!/usr/bin/env Rscript
# Bimodality / sign-label-ambiguity diagnostic for Omega_B[1,2]
# Usage: Rscript scripts/diagnostic_bimodality.R <out_dir>
#   out_dir: run output directory containing SUMMARY.txt and one_fit_n<N>_ci.rds
#   If omitted, falls back to hardcoded n=5 path for local development.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) {
  out_dir <- args[1L]
} else {
  out_dir <- "outputs/ci/phase0_n5_run26435798026"
}

summary_path <- file.path(out_dir, "SUMMARY.txt")
if (!file.exists(summary_path)) stop("Missing SUMMARY.txt: ", summary_path)

lines <- readLines(summary_path)
n_line <- grep("^N:", lines, value = TRUE)
if (length(n_line) == 0L) stop("N: line not found in SUMMARY.txt")
N <- as.integer(sub("^N:\\s*", "", n_line[1L]))
if (is.na(N)) stop("Could not parse N from SUMMARY.txt")

rds_path <- file.path(out_dir, sprintf("one_fit_n%d_ci.rds", N))
if (!file.exists(rds_path)) stop("Missing RDS: ", rds_path)

bundle   <- readRDS(rds_path)
rho_all  <- bundle$rho_B_posterior
n_chains <- bundle$fit_settings$chains
n_iter   <- bundle$fit_settings$sampling

if (is.null(rho_all)) stop("rho_B_posterior is NULL (fit may have crashed)")
if (length(rho_all) != n_chains * n_iter)
  stop(sprintf("Expected %d draws, got %d", n_chains * n_iter, length(rho_all)))

chain_draws <- lapply(seq_len(n_chains),
  function(ch) rho_all[((ch-1)*n_iter+1):(ch*n_iter)])
cq <- function(x, p) quantile(x, p, names=FALSE)
cstats <- lapply(chain_draws, function(x) c(med=median(x), lo=cq(x,.025), hi=cq(x,.975)))

png_path <- file.path(out_dir, "diagnostic_bimodality_pairs.png")
png(png_path, width=900, height=420, res=120)
par(mfrow=c(1L, n_chains+1L), mar=c(4,4,3,1))
plot(density(rho_all), main="Omega_B[1,2] all chains",
     xlab=expression(rho[B]), lwd=2); abline(v=0, lty=2, col="grey60")
cols <- c("steelblue","tomato","forestgreen","goldenrod")[seq_len(n_chains)]
for (ch in seq_len(n_chains)) {
  plot(density(chain_draws[[ch]]),
       main=sprintf("Chain %d  med=%+.3f", ch, cstats[[ch]]["med"]),
       xlab=expression(rho[B]), col=cols[ch], lwd=2)
  abline(v=0, lty=2, col="grey60")
}
dev.off()
cat("Pairs plot:", png_path, "\n")

signs     <- sapply(cstats, function(s) sign(s["med"]))
med_range <- diff(range(sapply(cstats, function(s) s["med"])))
rhat      <- bundle$omega_B_summary$rhat
ess       <- bundle$omega_B_summary$ess_bulk
verdict   <- if (length(unique(signs)) > 1L) "STRONGLY BIMODAL" else
  if (med_range > .40 || (!is.null(rhat) && rhat > 1.20)) "WEAKLY BIMODAL" else
  if (!is.null(rhat) && rhat < 1.10 && med_range < .20) "UNIMODAL" else
  "INSUFFICIENT EVIDENCE"

chain_lines <- sapply(seq_len(n_chains), function(ch) {
  s <- cstats[[ch]]
  sprintf("  Chain %d: median=%+.3f  95%% CrI [%+.3f, %+.3f]", ch, s["med"], s["lo"], s["hi"])
})

txt <- c(
  "PER-CHAIN MEDIANS:",
  chain_lines,
  sprintf("  (overall: Rhat=%.3f  ESS_bulk=%.0f)", rhat, ess),
  "",
  "PAIRS PLOT INTERPRETATION:",
  "  Bundle stores only Omega_B[1,2] draws; M[2,k] / tau_B not available.",
  "  Hypothesis TRUE => opposite-sign chain medians; bimodal overall density",
  "  (two peaks straddling 0); per-chain density unimodal but at opposite modes.",
  "  At n=48: M[2,2] (log-boost, biomarker-2) anti-correlated with Omega_B[1,2]",
  "  in joint scatter is the defining feature of the sign-flip.",
  "",
  paste("VERDICT:", verdict),
  "",
  "CAVEAT:",
  "  n=5 is weakly informative; bimodality may show only partial chain separation",
  "  (Rhat 1.1-1.5). For n=48 look for: (i) bimodal Omega_B[1,2] marginal;",
  "  (ii) M[2,2] anti-correlated with Omega_B[1,2]; (iii) Rhat > 1.1 with",
  "  ESS_bulk < 200 despite adequate iteration count."
)

verdict_path <- file.path(out_dir, "BIMODALITY_VERDICT.txt")
writeLines(txt, verdict_path)
cat("Verdict:", verdict_path, "\n")
cat("VERDICT:", verdict, "\n")
