#' @title Module-based Dynamic Network Biomarker (MDNB) Analysis
#' @description
#' Performing module-based dynamic network biomarker analysis using PPI network
#' @author Zaoqu Liu, Chuhan Zhang; Email: liuzaoqu@163.com
#'
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param state A dataframe with two columns: sample names and group information.
#' @param state.levels A vector for state sequence (e.g., cell groups, time points, treatment conditions).
#' @param cor.method Specifies the method for correlation analysis.
#' @param ppi Protein-protein interaction network; background network.
#' @param min.combined.score Minimum combined score for determining protein-protein interaction.
#' @param PCC.min Remove genes if their highest correlation with any other gene is less than this value.
#' @param PCC.module A gene must maintain at least this minimum correlation with the seed gene
#'   in every group to be included in the module.
#' @param numeber.module.QI Number of gene modules selected for QI calculation.
#'
#' @return A list containing:
#' \itemize{
#'   \item DNB.score: The DNB score values (QI values) calculated for different groups.
#'   \item DNB.genes: Genes from the module with highest CI in the critical state.
#'   \item CI_all: Module contribution index (CI) details for each state.
#'   \item Gene_module: Generated co-expression gene modules.
#'   \item QIplot: Line chart showing QI changes across states.
#'   \item metadata: Additional information about the analysis.
#' }
#'
#' @references
#' Li L, Xu Y, Yan L, et al. Dynamic network biomarker factors orchestrate cell-fate
#' determination at tipping points during hESC differentiation. Innovation (Camb).
#' 2022 Dec 20;4(1):100364. doi: 10.1016/j.xinn.2022.100364
#'
#' @export
MDNB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    ppi = ppi_h,
    min.combined.score = 900,
    PCC.min = 0.02,
    PCC.module = 0.02,
    numeber.module.QI = 10) {
  cat("Publish details: Li L et al. (2022) Innovation. DOI: 10.1016/j.xinn.2022.100364\n")
  cat("\n")

  cat("+++ Performing module-based dynamic network biomarker analysis...\n")
  cat("\n")

  # ============================================================================
  # 1. DATA PREPARATION
  # ============================================================================

  if (ncol(state) != 2) {
    stop("Number of state columns must be two: sample names and group information!")
  }

  colnames(state) <- c("ID", "State")
  expr <- expr[, colnames(expr) %in% state$ID]
  state <- state[match(colnames(expr), state$ID), ]

  # Create state index
  tmp <- aggregate(
    state$ID,
    by = list(group = state$State),
    FUN = function(x) paste(x, collapse = ",")
  )
  state_idx <- as.matrix(tmp)
  state_idx <- state_idx[order(factor(state_idx[, 1], levels = state.levels)), ]

  # Prepare PPI network
  # Use base R aggregate instead of data.table syntax to avoid dependency issues
  ppim <- ppi[, 1:2]
  ppim <- aggregate(G1 ~ G2, data = ppim, FUN = function(x) paste(x, collapse = ", "))
  colnames(ppim)[2] <- "G1_combined"

  # ============================================================================
  # 2. NORMALIZE EXPRESSION AND COMPUTE SD FOR EACH STATE
  # ============================================================================

  cat("+++ Normalizing expression and computing SD for each state...\n")

  for (d in 1:nrow(state_idx)) {
    long_str <- state_idx[d, 2]
    state_cols <- unlist(strsplit(long_str, split = ","))
    scaled_data <- t(scale(t(expr[, state_cols])))
    scaled_data[is.nan(scaled_data)] <- 0
    expr[, state_cols] <- scaled_data
  }
  expr <- as.data.frame(expr)

  # Compute standard deviations
  ssd <- matrix(0, nrow = nrow(expr), ncol = nrow(state_idx))
  rownames(ssd) <- rownames(expr)
  colnames(ssd) <- state_idx[, 1]

  for (i in 1:nrow(state_idx)) {
    long_str <- state_idx[i, 2]
    state_cols <- unlist(strsplit(long_str, split = ","))
    ssd[, i] <- apply(expr[, state_cols], 1, sd, na.rm = TRUE)
  }

  # ============================================================================
  # 3. SCREEN GENES IN PPI
  # ============================================================================

  cat("+++ Screening genes in PPI network...\n")
  feature <- rownames(expr)
  filtered_genes <- feature[feature %in% ppim[, 1]]

  cat("Original gene number:", length(feature), "→ Genes in PPI:", length(filtered_genes), "\n")

  # ============================================================================
  # 4. COMPUTE PAIRWISE PCC FOR ALL GENE PAIRS
  # ============================================================================

  cat("+++ Computing pairwise PCC for all gene pairs within groups...\n")
  
  # Check if C++ acceleration is available (REQUIRED!)
  has_cpp <- exists("fast_cor_cpp", where = asNamespace("EasyDNB"), mode = "function")
  
  if (!has_cpp) {
    stop("C++ acceleration not available! Please reinstall EasyDNB with:\n",
         "  cd /Users/liuzaoqu/Desktop/develop/EasyDNB\n",
         "  R CMD INSTALL .",
         call. = FALSE)
  }
  
  cat("    ✓ C++ acceleration ACTIVE (10-20x faster)\n")
  
  allpcc <- list()
  total_groups <- nrow(state_idx)
  pb <- txtProgressBar(min = 0, max = total_groups, style = 3, width = 50)

  for (t in 1:total_groups) {
    long_str <- state_idx[t, 2]
    state_cols <- unlist(strsplit(long_str, split = ","))
    expr_m <- expr[filtered_genes, state_cols]

    if (has_cpp) {
      # Use C++ acceleration
      pcc_matrix <- fast_cor_cpp(as.matrix(expr_m), method = cor.method)
      pcc_matrix <- abs(pcc_matrix)
    } else {
      # Fallback to R version
      pcc_matrix <- WGCNA::cor(t(expr_m), use = "p", method = cor.method)
      pcc_matrix[is.na(pcc_matrix)] <- 0
      pcc_matrix <- abs(pcc_matrix)
    }

    rownames(pcc_matrix) <- filtered_genes
    colnames(pcc_matrix) <- filtered_genes

    allpcc[[t]] <- pcc_matrix
    names(allpcc)[t] <- paste0("G", t)

    setTxtProgressBar(pb, t)
  }
  close(pb)

  # ============================================================================
  # 5. GENERATE GENE MODULES
  # ============================================================================

  cat("\n+++ Generating gene modules...\n")

  # Find seed genes
  max_pcc <- apply(allpcc[[1]], 1, max)
  seed_genes <- names(max_pcc[max_pcc >= PCC.min])

  cat("Number of seed genes (max PCC >=", PCC.min, "):", length(seed_genes), "\n")

  # Build modules
  modules <- list()
  pb <- txtProgressBar(min = 0, max = length(seed_genes), style = 3, width = 50)

  for (i in seq_along(seed_genes)) {
    seed_gene <- seed_genes[i]

    # Get PPI neighbors
    ppi_str <- ppim$G1_combined[ppim$G2 == seed_gene]
    if (length(ppi_str) > 0) {
      ppi_neighbors <- unlist(strsplit(ppi_str, split = ", "))
      ppi_neighbors <- intersect(ppi_neighbors, filtered_genes)

      if (length(ppi_neighbors) > 0) {
        # Check if neighbors maintain correlation across all groups
        valid_neighbors <- sapply(ppi_neighbors, function(neighbor) {
          all(sapply(allpcc, function(pcc_mat) {
            pcc_mat[seed_gene, neighbor] >= PCC.module
          }))
        })

        valid_neighbors <- ppi_neighbors[valid_neighbors]

        if (length(valid_neighbors) > 0) {
          modules[[i]] <- list(
            name = seed_gene,
            members = c(seed_gene, valid_neighbors)
          )
        }
      }
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Remove NULL modules
  modules <- modules[!sapply(modules, is.null)]

  cat("\nNumber of modules constructed:", length(modules), "\n")

  # ============================================================================
  # 6. CALCULATE CONTRIBUTION INDEX (CI) FOR EACH MODULE
  # ============================================================================

  cat("+++ Calculating Contribution Index (CI) for each module...\n")

  CI <- list()

  for (t in seq_along(allpcc)) {
    cat(sprintf("Group %d/%d: %s\n", t, length(allpcc), names(allpcc)[t]))

    cor_matrix <- allpcc[[t]]
    CI[[t]] <- matrix(0, nrow = length(modules), ncol = 6)

    pb_module <- txtProgressBar(min = 0, max = length(modules), style = 3, width = 50)

    for (i in seq_along(modules)) {
      module_gene_names <- modules[[i]]$members

      # 1. Module number and size
      CI[[t]][i, 1] <- i
      CI[[t]][i, 2] <- length(module_gene_names)

      # 2. Average SD
      CI[[t]][i, 4] <- mean(ssd[module_gene_names, t], na.rm = TRUE)

      # 3. Internal correlation (R_in)
      if (length(module_gene_names) > 1) {
        internal_genes <- module_gene_names[-1] # Exclude seed gene
        cor_submatrix <- cor_matrix[internal_genes, internal_genes, drop = FALSE]
        CI[[t]][i, 5] <- mean(cor_submatrix[lower.tri(cor_submatrix)], na.rm = TRUE)
      }

      # 4. External correlation (R_out)
      module_external_genes <- setdiff(filtered_genes, module_gene_names)

      if (length(modules[[i]]$members) > 1) {
        external_genes <- module_external_genes

        if (length(external_genes) > 0) {
          internal_gene_names <- module_gene_names[-1]
          cor_submatrix <- cor_matrix[internal_gene_names, external_genes, drop = FALSE]
          CI[[t]][i, 6] <- mean(cor_submatrix, na.rm = TRUE)
        }
      }

      # 5. Calculate CI
      module_size <- length(modules[[i]]$members)
      r_out <- CI[[t]][i, 6]

      # Handle invalid R_out to prevent division by zero or NA
      if (is.na(r_out) || is.nan(r_out) || r_out == 0 || abs(r_out) < 1e-10) {
        r_out <- 1e-10 # Use small value instead of 0
      }

      CI[[t]][i, 3] <- as.numeric(sqrt(module_size) * as.numeric(CI[[t]][i, 4]) *
        as.numeric(CI[[t]][i, 5]) / r_out)

      setTxtProgressBar(pb_module, i)
    }
    close(pb_module)

    colnames(CI[[t]]) <- c("Module", "Size", "CI", "SD_ave", "R_in", "R_out")
    rownames(CI[[t]]) <- paste0("Module", 1:length(modules))

    CI[[t]][is.na(CI[[t]])] <- 0
    CI[[t]][is.infinite(CI[[t]])] <- 0
    CI[[t]] <- as.data.frame(CI[[t]])
  }

  names(CI) <- state.levels
  names(modules) <- paste0("Module", 1:length(modules))

  # Create module names matrix
  module_names_matrix <- matrix(
    sapply(modules, function(x) x[["name"]]),
    ncol = 1
  )
  rownames(module_names_matrix) <- paste0("Module", 1:length(modules))
  colnames(module_names_matrix) <- "Module_Core_Gene"

  # ============================================================================
  # 7. CALCULATE QUANTITATIVE INDEX (QI)
  # ============================================================================

  cat("+++ Calculating Quantitative Index (QI)...\n")
  ng <- nrow(state_idx)
  topmCI <- vector("list", ng)
  avetopCI <- numeric(ng)

  for (t in seq_len(ng)) {
    ci_vec <- as.numeric(CI[[t]][, 3])
    sorted <- sort(ci_vec, decreasing = TRUE)
    n_top <- min(numeber.module.QI, length(ci_vec))
    avetopCI[t] <- mean(sorted[1:n_top])
    topmCI[[t]] <- sorted[1:n_top]
  }

  names(topmCI) <- state.levels
  QI <- avetopCI
  QI <- as.matrix(avetopCI)
  rownames(QI) <- state.levels
  colnames(QI) <- "QI"

  # ============================================================================
  # 8. DRAW QI CHANGE CHART
  # ============================================================================

  cat("+++ Drawing QI change chart...\n")
  QI_df <- data.frame(
    State = rownames(QI),
    DNB.score = QI[, 1]
  )
  QI_df$State <- factor(QI_df$State, levels = state.levels)

  p <- ggplot2::ggplot(QI_df, ggplot2::aes(x = State, y = DNB.score, group = 1)) +
    ggplot2::geom_line(color = "#41A98E", linewidth = 2) +
    ggplot2::geom_point(shape = 16, size = 4, color = "#ED6355") +
    ggplot2::labs(
      x = "Groups",
      y = "DNB Score (QI)",
      title = "Module-based DNB Scores"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 1),
      legend.position = "none",
      axis.text = ggplot2::element_text(size = 12, colour = "black"),
      axis.title = ggplot2::element_text(face = "bold", colour = "black", size = 14),
      plot.title = ggplot2::element_text(
        face = "bold", colour = "black",
        size = 14, hjust = 0.5
      )
    )

  print(p)

  cat("\n")

  # ============================================================================
  # 9. CALCULATE DNB GENES
  # ============================================================================

  cat("+++ Calculating DNB genes...\n")

  # Find the state with largest DNB score difference compared to adjacent states
  score_diff_prev <- c(0, diff(QI_df$DNB.score))
  score_diff_next <- c(-diff(QI_df$DNB.score), 0)

  peak_excess_scores <- ifelse(score_diff_prev > 0 & score_diff_next > 0,
    score_diff_prev + score_diff_next, 0
  )

  names(peak_excess_scores) <- QI_df$State

  peak_scores_exclude_first <- peak_excess_scores[-1]
  key_state <- names(which.max(peak_scores_exclude_first))

  # Get module with highest CI in critical state
  ci_matrix <- CI[[key_state]]
  ci_matrix_sorted <- ci_matrix[order(ci_matrix$CI, decreasing = TRUE), ]

  top_module <- ci_matrix_sorted$Module[1]
  DNB_genes <- modules[[top_module]][["members"]]

  cat(paste0("+++ ", key_state, " is the critical state...\n"))
  cat(paste0("+++ ", length(DNB_genes), " genes (DNB module) involved...\n"))
  cat("\n")
  cat("+++ Done!\n")

  # ============================================================================
  # 10. RETURN RESULTS
  # ============================================================================

  return(list(
    DNB.score = QI_df,
    DNB.genes = DNB_genes,
    CI_all = CI,
    Gene_module = modules,
    QIplot = p,
    metadata = list(
      n_modules = length(modules),
      critical_state = key_state,
      top_module = top_module,
      method = "MDNB"
    )
  ))
}
