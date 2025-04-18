# Example: Process Shigella data for use with serodynamics models

# Load required libraries
library(readxl)
library(dplyr)
library(shigella)  # Assuming your package is installed
library(serodynamics)

# Set seed for reproducibility
set.seed(123)

# Load example data from the package's extdata folder
file_path <- system.file("extdata", "3.8.2024 Compiled Shigella datav2.xlsx", package = "shigella")

df <- readxl::read_excel(file_path, sheet = "Compiled")

# Process IgA anti-IpaB antibody data from SOSAR study
dL_clean_ipab <- process_shigella_data(
  data         = df,
  study_filter = "SOSAR",
  antigen      = n_ipab_MFI
)

# View first few rows of cleaned case data
print(head(dL_clean_ipab))


