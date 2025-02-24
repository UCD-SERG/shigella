# create function for generating specific region data
create_xs_data <- function(df, filter_countries, filter_antigen_iso, value_col) {
  df %>%
    # First, create/rename columns
    mutate(
      id = sid,
      Country = site_name,
      study = study_name,
      age = age,
      antigen_iso = factor(isotype_name),
      value = {{ value_col }
      } # value_col is specified by the user, e.g. n_ipab_MFI
    ) %>%
    # Then filter by the desired catchment(s) and antigen iso(s)
    filter(
      Country %in% filter_countries,
      antigen_iso %in% filter_antigen_iso
    ) %>%
    # Create age categories
    mutate(
      ageCat = factor(case_when(
        age < 5 ~ "<5",
        age >= 5 & age <= 15 ~ "5-15",
        age > 15 ~ "16+"
      ))
    ) %>%
    # Optionally select only the needed columns
    select(id, Country, study, age, antigen_iso, value, ageCat)
}