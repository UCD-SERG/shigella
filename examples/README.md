# Examples Directory

This directory contains example scripts demonstrating how to use the `shigella` package functions for various analysis workflows.

## Available Examples

### 1. Data Processing Example (`01_data_processing_example.R`)
Demonstrates how to:
- Load and structure Shigella study data
- Process longitudinal data for specific studies and antigens
- Create cross-sectional datasets for geographic regions
- Prepare data for serocalculator analysis
- Generate summary statistics by age groups

**Usage:**
```r
source("examples/01_data_processing_example.R")
```

### 2. Simulation Workflow Example (`02_simulation_workflow_example.R`)
Demonstrates how to:
- Set up antibody decay curve parameters
- Define biological and measurement noise parameters
- Run seroincidence simulations
- Extract and summarize simulation results
- Evaluate coverage probabilities

**Usage:**
```r
source("examples/02_simulation_workflow_example.R")
```

### 3. Complete Analysis Workflow (`03_complete_analysis_workflow.R`)
Demonstrates a complete end-to-end analysis including:
- Multi-region data processing
- Regional comparisons and visualizations
- Age distribution analysis
- Antibody level distributions
- Incidence rate estimation and visualization
- Sample size and power analysis

**Usage:**
```r
source("examples/03_complete_analysis_workflow.R")
```

## Notes

- All examples use **mock/placeholder data** since the real Shigella dataset cannot be included in the repository due to privacy restrictions.
- When real data becomes available in the `data/` directory, simply replace the mock data creation sections with actual data loading.
- The examples are designed to work seamlessly once real `.rda` files are added to the `data/` directory.

## Requirements

Make sure all required packages are installed:
```r
install.packages(c("dplyr", "tibble", "ggplot2", "tidyr"))

# Install from GitHub
remotes::install_github("UCD-SERG/serodynamics")
remotes::install_github("UCD-SERG/serocalculator")
```

## Running Examples

You can run examples in several ways:

1. **Source the entire script:**
   ```r
   source("examples/01_data_processing_example.R")
   ```

2. **Run interactively:**
   Open the R file in RStudio and run sections step-by-step using Ctrl+Enter (Cmd+Enter on Mac).

3. **Run from command line:**
   ```bash
   Rscript examples/01_data_processing_example.R
   ```

## Customization

Feel free to modify these examples for your specific analysis needs. Each script is extensively commented to help you understand what each section does.
