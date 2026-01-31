#' Prepare Data Frame for Serocalculator Analysis
#'
#' Standardizes column names and sets required attributes for compatibility
#' with the `serocalculator` package, which requires specific column naming
#' conventions for age and antibody measurements.
#'
#' @param df A data frame containing cross-sectional antibody data with at
#'   minimum an age column and an antibody measurement column.
#' @param age_col Character string specifying the name of the age column in
#'   the input data frame. Default is "age".
#' @param value_col Character string specifying the name of the antibody
#'   measurement column. Default is "value".
#'
#' @return The input data frame with standardized column names and serocalculator
#'   attributes set. The age column is renamed to "age" if necessary, and
#'   attributes "age_var" and "value_var" are added for serocalculator recognition.
#'
#' @details
#' This function ensures that data frames are properly formatted for use with
#' the `serocalculator` package. It verifies that serocalculator can detect the
#' value column by checking internal package functions and provides appropriate
#' warnings or messages.
#'
#' @examples
#' \dontrun{
#' # Create mock antibody data
#' mock_data <- data.frame(
#'   id = 1:100,
#'   age = runif(100, 0.5, 60),
#'   antibody_level = rlnorm(100, log(1000), 1)
#' )
#'
#' # Prepare for serocalculator
#' prepared_data <- prepare_df_for_serocalculator(
#'   df = mock_data,
#'   age_col = "age",
#'   value_col = "antibody_level"
#' )
#' }
#'
#' @importFrom dplyr rename all_of %>%
#' @export
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