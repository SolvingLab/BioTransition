#' @title Sample-specific perturbation network based on PPI network
#' @description
#' Performing sample-specific perturbation network analysis based on PPI network.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param ref.samples Samples for constructing the background reference network; Column ID or names.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param ppi Protein-protein interaction network; background network.
#' @param min.combined.score Minimum combined score for determining protein-protein interaction.
#' @param nCores The number of cores will be used.
#' @import future
#' @import magrittr
#' @export
SSPN1 <- function(
    expr,
    ref.samples,
    cor.method = "pearson",
    p.adjust.method = "BH",
    ppi = ppi_h,
    min.combined.score = 900,
    nCores = parallel::detectCores() - 10) {
  cat("Publish details: Liu X et al. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016 Dec 15;44(22):e164. doi: 10.1093/nar/gkw772\n")
  cat("\n")

  cat("+++ Performing sample-specific perturbation network analysis...\n")
  cat("\n")
  
  # Auto-detect C++ availability (silently)
  use_cpp <- cor.method == "pearson" && 
             exists("fast_sspn_batch") && 
             is.function(fast_sspn_batch)
  if (is.numeric(ref.samples)) {
    case.samples <- colnames(expr)[-ref.samples]
  } else {
    case.samples <- colnames(expr)[!colnames(expr) %in% ref.samples]
  }

  ppi2 <- ppi[ppi$G1 %in% rownames(expr) & ppi$G2 %in% rownames(expr) & ppi$combined_score >= min.combined.score, -3]
  expr <- expr[unique(c(ppi2$G1)), ]
  cat(paste0("+++ ", nrow(expr), " intersected genes between expression matrix and PPI network...\n"))
  cat("\n")

  ppi3 <- as.data.frame(t(apply(ppi2, 1, sort))) %>% dplyr::distinct(V1, V2, .keep_all = T)
  cat(paste0("+++ ", nrow(ppi3), " edges in PPI network...\n"))
  cat("\n")

  cat("+++ Calculating background network using reference samples...\n")
  cat("\n")
  
  refExpr <- expr[, ref.samples]
  
  if (use_cpp) {
    ref_cor_result <- fast_cor_pval_cpp(t(as.matrix(refExpr)), method = "pearson")
    refCor <- ref_cor_result$r
    rownames(refCor) <- colnames(refCor) <- rownames(expr)
  } else {
    refCor <- suppressWarnings(CorandPval(t(refExpr),
      method = cor.method,
      p.adjust.method = p.adjust.method
    ))$r
  }
  
  ppi3$refCor <- sapply(1:nrow(ppi3), function(i) refCor[ppi3$V1[i], ppi3$V2[i]])

  cat("+++ Calculating sample-specific perturbation network for each case sample...\n")
  cat("\n")
  
  if (use_cpp) {
    # C++ batch processing
    
    # Convert to 0-based indices
    case_indices_0based <- match(case.samples, colnames(expr)) - 1
    ref_indices_0based <- match(ref.samples, colnames(expr)) - 1
    
    # Map PPI genes to indices
    gene_to_idx <- setNames(0:(nrow(expr) - 1), rownames(expr))
    ppi_idx1 <- gene_to_idx[as.character(ppi3$V1)]
    ppi_idx2 <- gene_to_idx[as.character(ppi3$V2)]
    
    # Batch calculation
    cpp_results <- fast_sspn_batch(
      expr_mat = as.matrix(expr),
      ref_indices = as.integer(ref_indices_0based),
      case_indices = as.integer(case_indices_0based),
      ppi_gene1 = as.character(ppi3$V1),
      ppi_gene2 = as.character(ppi3$V2),
      gene_names = rownames(expr)
    )
    
    # Adjust p-values for each sample
    caseSSPN.list <- lapply(1:length(case.samples), function(i) {
      df <- as.data.frame(cpp_results[[i]])
      
      # Adjust p-values
      if (p.adjust.method == "BH") {
        df$caseFDR <- fast_bh_adjust(df$caseP)
      } else {
        df$caseFDR <- p.adjust(df$caseP, method = p.adjust.method)
      }
      
      # Add gene names
      df$V1 <- ppi3$V1
      df$V2 <- ppi3$V2
      
      return(df[, c("V1", "V2", "corPert", "caseZ", "caseP", "caseFDR")])
    })
    names(caseSSPN.list) <- case.samples
    
  } else {
    # R version with parallel processing
    tmp_fun <- function(id) {
      caseExpr <- cbind(refExpr, case = expr[, id])
      caseCor <- suppressWarnings(CorandPval(t(caseExpr), method = cor.method, p.adjust.method = p.adjust.method))$r
      caseCor2 <- purrr::map_vec(1:nrow(ppi3), ~ caseCor[ppi3$V1[.x], ppi3$V2[.x]])
      casePert <- caseCor2 - ppi3$refCor
      caseZ <- casePert / ((1 - ppi3$refCor^2) / (length(ref.samples) - 1))
      caseP <- purrr::map_vec(caseZ, function(x) {
        2 * stats::pnorm(-abs(x))
      })
      caseP[is.na(caseP) | is.nan(caseP)] <- 1
      caseFDR <- stats::p.adjust(caseP, p.adjust.method)
      return(data.frame(V1 = ppi3$V1, V2 = ppi3$V2, corPert = casePert, caseZ = caseZ, caseP = caseP, caseFDR = caseFDR))
    }
    
    plan(multisession, workers = nCores)
    caseSSPN.list <- furrr::future_map(case.samples, ~ tmp_fun(.x), .progress = TRUE)
    names(caseSSPN.list) <- case.samples
    future::plan(future::sequential)
  }

  case.SSPN.Pert.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "corPert"]))
  case.SSPN.P.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "caseP"]))
  case.SSPN.FDR.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "caseFDR"]))
  rownames(case.SSPN.Pert.matrix) <- rownames(case.SSPN.P.matrix) <-
    rownames(case.SSPN.FDR.matrix) <- paste0(ppi3$V1, "_", ppi3$V2)

  cat("\n")
  cat("+++ Done!\n")
  future::plan(future::sequential)

  return(list(
    case.SSPN.list = caseSSPN.list,
    refNet = ppi3,
    case.SSPN.Pert.matrix = case.SSPN.Pert.matrix,
    case.SSPN.P.matrix = case.SSPN.P.matrix,
    case.SSPN.FDR.matrix = case.SSPN.FDR.matrix,
    PPI.used = ppi3[, -3]
  ))
}

#' @title Sample-specific perturbation network based on customized network
#' @description
#' Performing sample-specific perturbation network analysis based on customized network.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param ref.samples Samples for constructing the background reference network; Column ID or names.
#' @param net A customized network with two columns.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param nCores The number of cores will be used.
#' @import future
#' @import magrittr
#' @export
SSPN2 <- function(
    expr,
    ref.samples,
    net,
    cor.method = "pearson",
    p.adjust.method = "BH",
    nCores = parallel::detectCores() - 10) {
  cat("Publish details: Liu X et al. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016 Dec 15;44(22):e164. doi: 10.1093/nar/gkw772\n")
  cat("\n")

  cat("+++ Performing sample-specific perturbation network analysis...\n")
  cat("\n")
  
  # Auto-detect C++ availability (silently)
  use_cpp <- cor.method == "pearson" && 
             exists("fast_sspn_batch") && 
             is.function(fast_sspn_batch)
  if (is.numeric(ref.samples)) {
    case.samples <- colnames(expr)[-ref.samples]
  } else {
    case.samples <- colnames(expr)[!colnames(expr) %in% ref.samples]
  }

  colnames(net) <- c("G1", "G2")
  net2 <- net[net$G1 %in% rownames(expr) & net$G2 %in% rownames(expr), ]
  expr <- expr[unique(c(net2$G1, net2$G2)), ]
  cat(paste0("+++ ", nrow(expr), " intersected genes between expression matrix and net network...\n"))
  cat("\n")

  net3 <- dplyr::distinct(net2, G1, G2, .keep_all = T)
  cat(paste0("+++ ", nrow(net3), " edges in net network...\n"))
  cat("\n")

  cat("+++ Calculating background network using reference samples...\n")
  cat("\n")

  refExpr <- expr[, ref.samples]
  
  if (use_cpp) {
    ref_cor_result <- fast_cor_pval_cpp(t(as.matrix(refExpr)), method = "pearson")
    refCor <- ref_cor_result$r
    rownames(refCor) <- colnames(refCor) <- rownames(expr)
  } else {
    refCor <- suppressWarnings(CorandPval(t(refExpr),
      method = cor.method,
      p.adjust.method = p.adjust.method
    ))$r
  }
  
  net3$refCor <- sapply(1:nrow(net3), function(i) refCor[net3$G1[i], net3$G2[i]])

  cat("+++ Calculating sample-specific perturbation network for each case sample...\n")
  cat("\n")
  
  if (use_cpp) {
    # C++ batch processing
    
    # Convert to 0-based indices
    case_indices_0based <- match(case.samples, colnames(expr)) - 1
    ref_indices_0based <- match(ref.samples, colnames(expr)) - 1
    
    # Map PPI genes to indices (0-based)
    gene_to_idx <- setNames(0:(nrow(expr) - 1), rownames(expr))
    ppi_idx1 <- gene_to_idx[as.character(net3$G1)]
    ppi_idx2 <- gene_to_idx[as.character(net3$G2)]
    
    # Batch calculation
    cpp_results <- fast_sspn_batch(
      expr_mat = as.matrix(expr),
      ref_indices = as.integer(ref_indices_0based),
      case_indices = as.integer(case_indices_0based),
      ppi_gene1 = as.character(net3$G1),
      ppi_gene2 = as.character(net3$G2),
      gene_names = rownames(expr)
    )
    
    # Adjust p-values for each sample
    caseSSPN.list <- lapply(1:length(case.samples), function(i) {
      df <- as.data.frame(cpp_results[[i]])
      
      # Adjust p-values
      if (p.adjust.method == "BH") {
        df$caseFDR <- fast_bh_adjust(df$caseP)
      } else {
        df$caseFDR <- p.adjust(df$caseP, method = p.adjust.method)
      }
      
      # Add gene names
      df$G1 <- net3$G1
      df$G2 <- net3$G2
      
      return(df[, c("G1", "G2", "corPert", "caseZ", "caseP", "caseFDR")])
    })
    names(caseSSPN.list) <- case.samples
    
  } else {
    tmp_fun <- function(id) {
      caseExpr <- cbind(refExpr, case = expr[, id])
      caseCor <- suppressWarnings(CorandPval(t(caseExpr), method = cor.method, p.adjust.method = p.adjust.method))$r
      caseCor2 <- purrr::map_vec(1:nrow(net3), ~ caseCor[net3$G1[.x], net3$G2[.x]])
      casePert <- caseCor2 - net3$refCor
      caseZ <- casePert / ((1 - net3$refCor^2) / (length(ref.samples) - 1))
      caseP <- purrr::map_vec(caseZ, function(x) {
        2 * stats::pnorm(-abs(x))
      })
      caseP[is.na(caseP) | is.nan(caseP)] <- 1
      caseFDR <- stats::p.adjust(caseP, p.adjust.method)
      return(data.frame(G1 = net3$G1, G2 = net3$G2, corPert = casePert, caseZ = caseZ, caseP = caseP, caseFDR = caseFDR))
    }
    
    plan(multisession, workers = nCores)
    caseSSPN.list <- furrr::future_map(case.samples, ~ tmp_fun(.x), .progress = TRUE)
    names(caseSSPN.list) <- case.samples
    future::plan(future::sequential)
  }

  case.SSPN.Pert.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "corPert"]))
  case.SSPN.P.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "caseP"]))
  case.SSPN.FDR.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "caseFDR"]))
  rownames(case.SSPN.Pert.matrix) <- rownames(case.SSPN.P.matrix) <-
    rownames(case.SSPN.FDR.matrix) <- paste0(net3$G1, "_", net3$G2)

  cat("\n")
  cat("+++ Done!\n")
  future::plan(future::sequential)

  return(list(
    case.SSPN.list = caseSSPN.list,
    refNet = net3,
    case.SSPN.Pert.matrix = case.SSPN.Pert.matrix,
    case.SSPN.P.matrix = case.SSPN.P.matrix,
    case.SSPN.FDR.matrix = case.SSPN.FDR.matrix,
    Net.used = net3[, -3]
  ))
}
