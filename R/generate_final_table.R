#' Generate Final Table from Simulation Results
#'
#' This function loops through a list of simulation results, extracts the required columns,
#' and combines them into a single data frame. It also adds the sample size for clarity.
#'
#' @param results_list A list of simulation results. Each element should contain an element named \code{est1}
#' which can be summarized.
#' @param sample_size An integer specifying the sample size used in the simulation.
#'
#' @return A data frame combining the selected columns (\code{incidence.rate}, \code{SE}, \code{CI.lwr}, \code{CI.upr})
#' from each simulation result along with the sample size and an index for tracking.
#' @examples
#' \dontrun{
#' # Suppose you have a list of simulation results
#' final_table <- generate_final_table(results_list, sample_size = 100)
#' }
#' @export
generate_final_table <- function(results_list, sample_size) {
  # Initialize an empty list to store the results
  summary_results <- list()
  
  # Loop through each result and extract the required columns
  for (i in 1:length(results_list)) {
    # Extract the summary for each result
    result_summary <- summary(results_list[[i]]$est1)
    
    # Select the required columns (ensure that result_summary is a data frame)
    extracted_columns <- result_summary %>%
      select(incidence.rate, SE, CI.lwr, CI.upr)
    
    # Add a column for the index (optional, for tracking)
    extracted_columns <- extracted_columns %>%
      mutate(index = i)
    
    # Append the extracted data to the list
    summary_results[[i]] <- extracted_columns
  }
  
  # Combine all results into a single data frame and add the sample size column for clarity
  final_table <- bind_rows(summary_results) %>%
    mutate(sample_size = sample_size)
  
  return(final_table)
}
