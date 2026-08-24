#!/usr/bin/env Rscript
# =====================================================================
# 01_make_xs_subset.R -- build the cross-sectional analysis set for Chapter 3
#
# Kristen's specification (meeting 2026-08-03):
#   - the population-based serosurvey arms are `highE` and `lowE`
#   - the `_hh` arms are household members enrolled alongside an index case
#     and must be left out
#   - the survey was longitudinal (baseline, 3 months, 6 months, later
#     follow-up); the Lancet Microbe analysis used the baseline visit, and so
#     do we
#
# WHAT THIS PRODUCES
#   sees_xs_subset.rds       long format, one row per person-biomarker
#   sees_xs_summary.txt      the filter cascade, for the methods section
#
# A NOTE ON WHICH VALUE COLUMN. `result` stores measurements below the limit
# of detection as 0, while `result.llod` replaces them with the limit itself.
# serocalculator's noise model integrates the response density from `y.low`
# upward, so a value of 0 is outside its support; `result.llod` is the column
# that matches the model. About 31% of HlyE IgG and 16% of HlyE IgA sit at the
# limit, so this is not a rare edge case.
#
#   Rscript R/01_make_xs_subset.R
#   ELISA=/path/to/elisa_clean.csv Rscript R/01_make_xs_subset.R
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

ELISA <- Sys.getenv("ELISA", "~/data/elisa_clean_2025-12-11.csv")
OUT <- Sys.getenv("OUT_DIR", "data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(..., "\n", sep = "")
hr <- function() cat(strrep("-", 70), "\n")

# ---- read only the columns we need ----------------------------------------
KEEP <- c(
  "pid", "index_id", "Country", "sees_site", "studyarm", "studyvisit",
  "antiso", "antigen", "isotype", "Age", "sampleage",
  "result", "result.llod", "llod", "sex", "catchment"
)

say("reading ", ELISA)
raw <- readr::read_csv(path.expand(ELISA),
  show_col_types = FALSE,
  col_select = any_of(KEEP), progress = FALSE
)
say("  ", format(nrow(raw), big.mark = ","), " rows, ", ncol(raw), " columns")

cascade <- tibble::tibble(
  step = character(), rows = integer(),
  people = integer(), note = character()
)
add <- function(d, step, note = "") {
  cascade <<- dplyr::bind_rows(cascade, tibble::tibble(
    step = step, rows = nrow(d), people = dplyr::n_distinct(d$pid), note = note
  ))
  d
}
raw <- add(raw, "all ELISA records")

# ---- 1. the two population-based survey arms ------------------------------
d <- raw |>
  filter(studyarm %in% c("highE", "lowE")) |>
  add(
    "population survey arms (highE, lowE)",
    "excludes _hh arms: household members of index cases"
  )

# ---- 2. baseline visit ------------------------------------------------------
d <- d |>
  filter(studyvisit == "baseline") |>
  add("baseline visit", "the survey was longitudinal; Lancet Microbe used baseline")

# ---- 3. the two HlyE isotypes ----------------------------------------------
d <- d |>
  filter(antiso %in% c("HlyE IgG", "HlyE IgA")) |>
  add("HlyE IgG and IgA", "the two biomarkers Chapter 2 models jointly")

# ---- 4. complete age and value ---------------------------------------------
d <- d |>
  filter(!is.na(Age), !is.na(result.llod), !is.na(llod)) |>
  add("complete age and measurement")

# ---- 5. both isotypes present ----------------------------------------------
# The joint likelihood pairs a person's two measurements, so a person with only
# one contributes nothing that the product form does not already give. Keeping
# them would make the two likelihoods answer slightly different questions.
both <- d |>
  count(pid, name = "n_iso") |>
  filter(n_iso == 2)
d <- d |>
  semi_join(both, by = "pid") |>
  add("both isotypes present", "required for the joint likelihood")

# ---- tidy ------------------------------------------------------------------
xs <- d |>
  transmute(
    id          = pid,
    index_id    = index_id,
    Country     = Country,
    site        = sees_site,
    arm         = studyarm,
    endemicity  = if_else(studyarm == "highE", "High", "Low"),
    age         = as.numeric(Age),
    antigen_iso = gsub(" ", "_", antiso), # "HlyE IgG" -> "HlyE_IgG"
    value       = as.numeric(result.llod),
    llod        = as.numeric(llod),
    sex         = sex
  ) |>
  arrange(arm, id, antigen_iso)

stopifnot(all(xs$antigen_iso %in% c("HlyE_IgG", "HlyE_IgA")))

# ---- report -----------------------------------------------------------------
hr()
say("FILTER CASCADE")
hr()
print(as.data.frame(cascade), row.names = FALSE)

hr()
say("ANALYSIS SET")
hr()
say(
  "  people ", format(n_distinct(xs$id), big.mark = ","),
  "   rows ", format(nrow(xs), big.mark = ",")
)
cat("\n")
print(xs |> count(arm, Country, name = "rows") |>
  tidyr::pivot_wider(names_from = Country, values_from = rows, values_fill = 0) |>
  as.data.frame(), row.names = FALSE)
cat("\n")
print(xs |> group_by(arm) |>
  summarise(
    people = n_distinct(id),
    age_min = min(age), age_q1 = quantile(age, .25),
    age_med = median(age), age_q3 = quantile(age, .75),
    age_max = max(age), .groups = "drop"
  ) |>
  as.data.frame(), row.names = FALSE)
cat("\n")
print(xs |> group_by(antigen_iso, arm) |>
  summarise(
    n = n(), at_llod = sum(value <= llod * 1.001),
    pct_at_llod = round(100 * mean(value <= llod * 1.001), 1),
    med = round(median(value), 3), .groups = "drop"
  ) |>
  as.data.frame(), row.names = FALSE)

# ---- the detection limits must match the noise parameters ------------------
hr()
say("LIMITS OF DETECTION BY COUNTRY")
say("  These must match y.low in the noise parameters, or the response density")
say("  will be integrated over the wrong support.")
hr()
print(xs |> distinct(Country, antigen_iso, llod) |>
  arrange(antigen_iso, Country) |> as.data.frame(), row.names = FALSE)

if (requireNamespace("serocalculator", quietly = TRUE)) {
  np <- serocalculator::example_noise_params_sees
  chk <- xs |>
    distinct(Country, antigen_iso, llod) |>
    left_join(np |> select(antigen_iso, Country, y.low),
      by = c("antigen_iso", "Country")
    ) |>
    mutate(match = abs(llod - y.low) < 1e-4)
  cat("\n")
  print(as.data.frame(chk), row.names = FALSE)
  if (all(chk$match, na.rm = TRUE) && !any(is.na(chk$y.low))) {
    say("\n  every limit matches example_noise_params_sees")
  } else {
    say("\n  *** a limit does not match. Do not proceed until this is resolved. ***")
  }
}

saveRDS(xs, file.path(OUT, "sees_xs_subset.rds"))
capture.output(
  {
    print(as.data.frame(cascade), row.names = FALSE)
    cat("\n")
    print(xs |> count(arm, Country) |> as.data.frame(), row.names = FALSE)
  },
  file = file.path(OUT, "sees_xs_summary.txt")
)

hr()
say("saved -> ", file.path(OUT, "sees_xs_subset.rds"))
say("        ", file.path(OUT, "sees_xs_summary.txt"))
