test_that("cDNB works with basic input", {
  set.seed(42)
  

  # Create minimal test data
  n_genes <- 100
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
  
  # Run cDNB (suppress output)
  result <- suppressMessages(
    cDNB(
      expr = expr,
      state = state,
      state.levels = c("A", "B", "C"),
      cor.method = "pearson",
      variation.method = "sd"
    )
  )
  
  # Check output structure

expect_type(result, "list")
  expect_true("DNB.score" %in% names(result))
  expect_true("DNB.genes" %in% names(result))
  
  # Check DNB.score structure
  expect_s3_class(result$DNB.score, "data.frame")
  expect_equal(nrow(result$DNB.score), 3)  # 3 states
  expect_true("CI" %in% colnames(result$DNB.score))
  
  # Check DNB.genes is character vector
  expect_type(result$DNB.genes, "character")
  expect_true(length(result$DNB.genes) > 0)
})

test_that("cDNB handles different correlation methods", {
  set.seed(42)
  
  n_genes <- 50
  n_samples <- 12
  
  expr <- matrix(
    rnorm(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(expr) <- paste0("Gene", 1:n_genes)
  colnames(expr) <- paste0("Sample", 1:n_samples)
  
  state <- data.frame(
    sample_id = colnames(expr),
    state = rep(c("X", "Y"), each = 6)
  )
  
  # Test with spearman correlation
  result <- suppressMessages(
    cDNB(
      expr = expr,
      state = state,
      state.levels = c("X", "Y"),
      cor.method = "spearman"
    )
  )
  
  expect_type(result, "list")
  expect_true("DNB.score" %in% names(result))
})

test_that("cDNB validates input correctly", {
  set.seed(42)
  
  expr <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:10)
  colnames(expr) <- paste0("Sample", 1:10)
  
  # Wrong number of columns in state
  state_wrong <- data.frame(
    sample_id = colnames(expr),
    state = rep(c("A", "B"), each = 5),
    extra = 1:10  # Extra column
  )
  
  expect_error(
    suppressMessages(
      cDNB(expr = expr, state = state_wrong, state.levels = c("A", "B"))
    )
  )
})
