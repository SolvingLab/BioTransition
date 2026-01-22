test_that("tDNB works with basic input", {
  set.seed(42)
  
  # Create minimal test data
  n_genes <- 80
  n_samples <- 15
  
  expr <- matrix(
    rnorm(n_genes * n_samples, mean = 10, sd = 2),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(expr) <- paste0("Gene", 1:n_genes)
  colnames(expr) <- paste0("Sample", 1:n_samples)
  
  state <- data.frame(
    sample_id = colnames(expr),
    state = rep(c("A", "B", "C"), each = 5)
  )
  
  # Run tDNB (suppress output)
  result <- suppressMessages(
    tDNB(
      expr = expr,
      state = state,
      state.levels = c("A", "B", "C"),
      cor.method = "pearson",
      variation.method = "sd",
      min.size = 5,
      max.size = 50
    )
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_true("DNB.score" %in% names(result))
  expect_true("DNB.genes" %in% names(result))
  expect_true("Tom" %in% names(result))  # tDNB specific
  expect_true("SFNet" %in% names(result))  # tDNB specific
  
  # Check DNB.score structure
  expect_s3_class(result$DNB.score, "data.frame")
  expect_equal(nrow(result$DNB.score), 3)  # 3 states
  expect_true("CI" %in% colnames(result$DNB.score))
  
  # Check DNB.genes is character vector
  expect_type(result$DNB.genes, "character")
  expect_true(length(result$DNB.genes) > 0)
})

test_that("tDNB returns TOM matrices for each state", {
  set.seed(123)
  
  n_genes <- 60
  n_samples <- 10
  
  expr <- matrix(
    rnorm(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(expr) <- paste0("Gene", 1:n_genes)
  colnames(expr) <- paste0("Sample", 1:n_samples)
  
  state <- data.frame(
    sample_id = colnames(expr),
    state = rep(c("X", "Y"), each = 5)
  )
  
  result <- suppressMessages(
    tDNB(
      expr = expr,
      state = state,
      state.levels = c("X", "Y"),
      min.size = 5,
      max.size = 40
    )
  )
  
  # Check TOM matrices
  expect_type(result$Tom, "list")
  expect_equal(length(result$Tom), 2)  # 2 states
  expect_true(all(c("X", "Y") %in% names(result$Tom)))
})
