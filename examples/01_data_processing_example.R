# Example 1: Data Processing Workflow
# This script demonstrates how to process raw Shigella study data
# using the shigella package functions

library(shigella)
library(dplyr)
library(tibble)

# -----------------------------------------------------------------------------
# 1. Create Mock Shigella Study Data
# -----------------------------------------------------------------------------
# Note: Real data will be loaded from data/ directory when available
# For now, we create placeholder data with the expected structure

set.seed(123)

mock_shigella_data <- tibble(
  study_name = rep(c("SOSAR", "GEMS"), each = 50),
  site_name = sample(c("MA USA", "Ghana", "Niger", "Sierra Leone"), 100, replace = TRUE),
  sid = rep(1:50, 2),
  isotype_name = sample(c("IgG", "IgA"), 100, replace = TRUE),
  timepoint = rep(c("V1", "V2"), 50),
  `Actual day` = rep(c(0, 14), 50) + rnorm(100, 0, 2),
  age = runif(100, 0.5, 60),
  ipab_MFI = rlnorm(100, log(1000), 1),
  sf2a_MFI = rlnorm(100, log(800), 1.2),
  sf3a_MFI = rlnorm(100, log(1200), 0.9),
  n_ipab_MFI = rlnorm(100, log(1000), 1),
  n_sf2a_MFI = rlnorm(100, log(800), 1.2),
  n_sf3a_MFI = rlnorm(100, log(1200), 0.9)
)

# -----------------------------------------------------------------------------
# 2. Process Longitudinal Data for a Specific Study and Antigen
# -----------------------------------------------------------------------------

processed_sosar <- process_shigella_data(
  data = mock_shigella_data,
  study_filter = "SOSAR",
  antigen = ipab_MFI
)

cat("\nProcessed SOSAR longitudinal data:\n")
print(head(processed_sosar, 10))

# -----------------------------------------------------------------------------
# 3. Create Cross-Sectional Data for Specific Regions
# -----------------------------------------------------------------------------

# Extract cross-sectional data for USA, IgG only
xs_usa_igg <- create_xs_data(
  df = mock_shigella_data,
  filter_countries = c("MA USA"),
  filter_antigen_iso = c("IgG"),
  value_col = n_ipab_MFI
)

cat("\nCross-sectional data for USA (IgG):\n")
print(head(xs_usa_igg, 10))

# Extract cross-sectional data for Ghana
xs_ghana_igg <- create_xs_data(
  df = mock_shigella_data,
  filter_countries = c("Ghana"),
  filter_antigen_iso = c("IgG"),
  value_col = n_ipab_MFI
)

cat("\nCross-sectional data for Ghana (IgG):\n")
print(head(xs_ghana_igg, 10))

# -----------------------------------------------------------------------------
# 4. Summary Statistics by Age Group
# -----------------------------------------------------------------------------

summary_by_age <- xs_usa_igg %>%
  group_by(ageCat) %>%
  summarise(
    n = n(),
    mean_value = mean(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

cat("\nSummary statistics by age category:\n")
print(summary_by_age)

# -----------------------------------------------------------------------------
# 5. Prepare Data for Serocalculator
# -----------------------------------------------------------------------------

xs_usa_prepared <- prepare_df_for_serocalculator(
  df = xs_usa_igg,
  age_col = "age",
  value_col = "value"
)

cat("\nData prepared for serocalculator analysis\n")
cat("Attributes:", paste(names(attributes(xs_usa_prepared)), collapse = ", "), "\n")
