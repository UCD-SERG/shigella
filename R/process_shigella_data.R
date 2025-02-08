process_shigella_data <- function(data, study_filter, antigen) {
  # Filter the data for the specific study
  filtered_data <- data %>%
    filter(study_name == study_filter)
  
  # Capture the column name of the antigen
  antigen_col <- ensym(antigen)
  
  # Manipulate and restructure the data
  processed_data <- filtered_data %>%
    select(isotype_name, sid, timepoint, `Actual day`, !!antigen_col) %>%
    mutate(
      index_id = sid,
      antigen_iso = isotype_name,
      visit = timepoint,
      timeindays = `Actual day`,
      result = !!antigen_col
    ) %>%
    group_by(index_id, antigen_iso) %>%
    arrange(visit) %>%
    mutate(visit_num = rank(visit, ties.method = "first")) %>%
    ungroup() %>%
    # Remove rows with NA in timeindays
    filter(!is.na(timeindays))
  
  return(processed_data)
}