# create table of incidence.rate of each region
create_incidence_table <- function(...) {
  # Capture input objects and their names
  est_list <- list(...)
  country_names <- names(est_list)
  
  # Extract incidence.rate from the summary() of each estimate object
  incidence_rates <- sapply(est_list, function(x) summary(x)$incidence.rate)
  
  # Create a tidy tibble
  incidence_table <- tibble(
    Country = country_names,
    Incidence_Rate = incidence_rates
  )
  
  return(incidence_table)
}