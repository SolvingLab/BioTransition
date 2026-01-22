test_that("PPI data is available and correctly formatted", {
  # Test human PPI
  data(ppi_h, package = "BioTransition")
  
  expect_s3_class(ppi_h, "data.frame")
  expect_true(ncol(ppi_h) >= 2)
  expect_true(nrow(ppi_h) > 0)
  
  # Test mouse PPI
  data(ppi_m, package = "BioTransition")
  
  expect_s3_class(ppi_m, "data.frame")
  expect_true(ncol(ppi_m) >= 2)
  expect_true(nrow(ppi_m) > 0)
})

test_that("PPI networks contain valid gene pairs", {
  data(ppi_h, package = "BioTransition")
  
  # Check no NA values in first two columns
  expect_false(any(is.na(ppi_h[, 1])))
  expect_false(any(is.na(ppi_h[, 2])))
  
  # Check gene names are character
  expect_type(ppi_h[[1]], "character")
  expect_type(ppi_h[[2]], "character")
})
