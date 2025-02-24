# Create function to clearly define cross-sectional data
prepare_df_for_serocalculator <- function(df, age_col = "age", value_col = "value") {
  # Ensure correct column names
  df <- df %>%
    rename(age = all_of(age_col))
  
  # Assign attributes for serocalculator
  attr(df, "age_var") <- "age"
  attr(df, "value_var") <- value_col
  
  # Check if serocalculator recognizes attributes
  get_value_var <- serocalculator:::get_value_var
  detected_value_var <- get_value_var(df)
  
  if (is.null(detected_value_var)) {
    warning("serocalculator did not detect the 'value' column. Check column naming.")
  } else {
    message("serocalculator recognized 'value' column: ", detected_value_var)
  }
  
  return(df)
}