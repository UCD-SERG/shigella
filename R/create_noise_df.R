create_noise_df <- function(country) {
  noise_df <- tibble(
    antigen_iso = c("IgG"),
    Country     = factor(c(country)), # Use input argument for Country
    y.low       = c(25),
    eps         = c(0.25),
    nu          = c(0.5),
    y.high      = c(200000)
  )
  
  return(noise_df)
}