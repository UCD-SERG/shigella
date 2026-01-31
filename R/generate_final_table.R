#' Generate Summary Table from Simulation Results
#'
#' Extracts and combines seroincidence estimates from multiple simulation
#' replicates into a single summary table for analysis and visualization.
#'
#' @param results_list A list of simulation results, where each element contains
#'   an `est1` component with seroincidence estimates (output from
#'   `simulate_seroincidence()`).
#' @param sample_size Integer. The sample size used in the simulations, added
#'   as a column for tracking and comparison across different sample sizes.
#'
#' @return A data frame with columns:
#'   \item{incidence.rate}{Estimated incidence rate (per person-year)}
#'   \item{SE}{Standard error of the estimate}
#'   \item{CI.lwr}{Lower bound of the 95\% confidence interval}
#'   \item{CI.upr}{Upper bound of the 95\% confidence interval}
#'   \item{index}{Simulation replicate number}
#'   \item{sample_size}{Sample size used in the simulation}
#'
#' @examples
#' \dontrun{
#' # Assuming simulation results exist
#' results <- simulate_seroincidence(
#'   dmcmc = mock_dmcmc,
#'   nrep = 200,
#'   n_sim = 100,
#'   observed = 0.1
#' )
#'
#' # Generate summary table
#' summary_table <- generate_final_table(
#'   results_list = results,
#'   sample_size = 200
#' )
#' }
#'
#' @importFrom dplyr select mutate bind_rows %>%
#' @export
# Define a function to generate final tables
generate_final_table <- function(results_list, sample_size) {
  # Initialize an empty list to store the results
  summary_results <- list()
  
  # Loop through each of the 100 results and extract the required columns
  for (i in 1:200) {
    # Extract the summary for each result
    result_summary <- summary(results_list[[i]]$est1)
    
    # Select the required columns
    extracted_columns <- result_summary %>%
      select(incidence.rate, SE, CI.lwr, CI.upr)
    
    # Add a column for the index (optional, for tracking)
    extracted_columns <- extracted_columns %>%
      mutate(index = i)
    
    # Append to the list
    summary_results[[i]] <- extracted_columns
  }
  
  # Combine all results into a single data frame
  final_table <- bind_rows(summary_results) %>%
    mutate(sample_size = sample_size) # Add sample size column for clarity
  
  return(final_table)
}
