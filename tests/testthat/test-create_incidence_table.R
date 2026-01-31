test_that("create_incidence_table extracts rates correctly", {
  # Create mock estimate objects
  mock_est1 <- list()
  class(mock_est1) <- "incidence.est"
  attr(mock_est1, "summary") <- data.frame(incidence.rate = 0.15)
  
  mock_est2 <- list()
  class(mock_est2) <- "incidence.est"
  attr(mock_est2, "summary") <- data.frame(incidence.rate = 0.22)
  
  # Mock summary method
  summary.incidence.est <- function(x) {
    attr(x, "summary")
  }
  
  # This test is conceptual - actual implementation depends on 
  # the structure of serocalculator/serodynamics estimate objects
  skip("Requires actual estimate objects from serocalculator")
})

test_that("create_incidence_table uses input names", {
  # Placeholder test for structure
  # In practice, would need real estimate objects
  skip("Requires actual estimate objects from serocalculator")
})

test_that("create_incidence_table returns tibble with correct columns", {
  # Placeholder test
  skip("Requires actual estimate objects from serocalculator")
})
