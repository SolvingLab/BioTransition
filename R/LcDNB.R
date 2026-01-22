#' @title Local Conventional DNB Analysis
#'
#' @description
#' Performs local conventional dynamic network biomarker (LcDNB) analysis using
#' protein-protein interaction (PPI) networks. This method extends the
#' conventional DNB approach by incorporating network topology information.
#'
#' @param expr A numeric matrix or data.frame with genes in rows and samples in
#'   columns. Row names should be gene symbols.
#' @param state A data.frame with exactly two columns: sample identifiers
#'   (matching column names in \code{expr}) and state/group labels.
#' @param state.levels A character vector specifying the order of states for
#'   analysis (e.g., time points or disease stages).
#' @param cor.method Character string specifying correlation method. One of
#'   \code{"pearson"} (default), \code{"spearman"}, or \code{"kendall"}.
#' @param p.adjust.method Character string specifying p-value adjustment
#'   method. One of \code{"BH"} (default), \code{"bonferroni"}, \code{"holm"},
#'   \code{"hochberg"}, \code{"hommel"}, \code{"BY"}, \code{"fdr"}, or
#'   \code{"none"}.
#' @param variation.method Character string specifying the method for
#'   calculating gene variation: \code{"sd"} (standard deviation, default) or
#'   \code{"cv"} (coefficient of variation).
#' @param min.first.neighbor.size Integer. Minimum number of first-order
#'   neighbors required for a gene to be considered as a center gene.
#'   Default: 3.
#' @param min.second.neighbor.size Integer. Minimum number of second-order
#'   neighbors required. Default: 1.
#' @param ppi A data.frame containing protein-protein interactions. Must have
#'   columns \code{G1}, \code{G2}, and \code{combined_score}. Default uses
#'   built-in human PPI network.
#' @param min.combined.score Numeric. Minimum STRING combined score for
#'   filtering PPI interactions. Default: 900 (high confidence).
#' @param percent Logical. If \code{TRUE} (default), use \code{top.p}
#'   percentage; if \code{FALSE}, use \code{top.n} absolute number.
#' @param top.n Integer. Number of top DNB genes when \code{percent = FALSE}.
#'   Default: 30.
#' @param top.p Numeric. Proportion of top DNB genes when \code{percent = TRUE}.
#'   Default: 0.05 (5\%).
#' @param AddModuleSize Logical. Whether to weight DNB score by module size.
#'   Default: \code{FALSE}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{DNB.score}{A data.frame with states and their DNB scores}
#'   \item{DNB.genes}{Character vector of identified DNB genes}
#'   \item{CI_all}{List of data.frames with CI scores for all center genes in
#'     each state}
#'   \item{Gene_module}{List of first-order neighbor genes for each center gene}
#'   \item{Candidate}{Data.frame summarizing candidate information per state}
#'   \item{Cor}{List of correlation matrices for each state}
#'   \item{V}{List of gene variation values for each state}
#'   \item{PPI.used}{Filtered PPI network used in analysis}
#'   \item{first.order.genes}{List of first-order neighbors}
#'   \item{second.order.genes}{List of second-order neighbors}
#' }
#'
#' @details
#' LcDNB analysis identifies critical transitions by combining gene expression
#' dynamics with PPI network topology. The algorithm:
#' \enumerate{
#'   \item Filters PPI network based on genes present in expression data
#'   \item Identifies center genes with sufficient network connections
#'   \item Calculates composite index (CI) using local network modules
#'   \item Identifies the critical state and DNB genes
#' }
#'
#' The composite index (CI) is calculated as:
#' CI = (SD_in * |PCC_in|) / |PCC_out|
#'
#' where SD_in is the mean standard deviation of module genes,
#' PCC_in is the mean correlation within the module, and
#' PCC_out is the mean correlation with external genes.
#'
#' @author Zaoqu Liu \email{liuzaoqu@@163.com}
#'
#' @seealso
#' \code{\link{cDNB}} for conventional DNB without PPI,
#' \code{\link{LDNB}} for landscape DNB,
#' \code{\link{ppi_h}} for human PPI network,
#' \code{\link{ppi_m}} for mouse PPI network
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' n_genes <- 200
#' n_samples <- 15
#'
#' expr <- matrix(
#'   rnorm(n_genes * n_samples, mean = 10, sd = 2),
#'   nrow = n_genes, ncol = n_samples
#' )
#' # Use gene names that exist in PPI network
#' data(ppi_h)
#' gene_names <- unique(c(ppi_h$G1, ppi_h$G2))[seq_len(n_genes)]
#' rownames(expr) <- gene_names
#' colnames(expr) <- paste0("Sample", seq_len(n_samples))
#'
#' state <- data.frame(
#'   sample_id = colnames(expr),
#'   state = rep(c("Normal", "Pre", "Disease"), each = 5)
#' )
#'
#' \donttest{
#' result <- LcDNB(
#'   expr = expr,
#'   state = state,
#'   state.levels = c("Normal", "Pre", "Disease"),
#'   min.combined.score = 400
#' )
#' result$DNB.score
#' head(result$DNB.genes)
#' }
#'
#' @export
LcDNB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    variation.method = "sd",
    min.first.neighbor.size = 3,
    min.second.neighbor.size = 1,
    ppi = ppi_h,
    min.combined.score = 900,
    percent = TRUE,
    top.n = 30,
    top.p = 0.05,
    AddModuleSize = FALSE) {

    ## Input validation
    .validate_expr(expr)
    .validate_state(state, colnames(expr))

    message("LcDNB: Local Conventional DNB Analysis")
    message("--------------------------------------")

    ## Prepare data
    colnames(state) <- c("ID", "state")
    state$state <- factor(state$state, levels = state.levels)
    state <- state[match(colnames(expr), state$ID), ]

    ## Split expression by state
    ddl <- lapply(state.levels, function(x) {
        expr[, state$ID[state$state == x], drop = FALSE]
    })
    names(ddl) <- state.levels

    ## Calculate correlations
    message("Calculating gene correlations...")
    corl <- lapply(ddl, function(d) {
        CorandPval(t(d), method = cor.method, p.adjust.method = p.adjust.method)
    })

    ## Calculate variations
    message("Calculating gene variations...")
    if (variation.method == "sd") {
        vl <- lapply(ddl, function(d) apply(d, 1, sd))
    } else if (variation.method == "cv") {
        vl <- lapply(ddl, function(d) {
            apply(d, 1, function(x) sd(x) / mean(x))
        })
    }

    data <- Map(function(x, y) list(r = x$r, v = y), corl, vl)

    ## Filter PPI network
    ppi2 <- ppi[ppi$G1 %in% rownames(expr) &
                ppi$G2 %in% rownames(expr) &
                ppi$combined_score >= min.combined.score, ]

    expr <- expr[unique(ppi2$G1), , drop = FALSE]
    message(sprintf("%d genes intersected with PPI network", nrow(expr)))

    ## Build local networks
    ppi3 <- ppi2[ppi2$G1 %in% names(which(table(ppi2$G1) >= min.first.neighbor.size)), ]
    ppil <- lapply(unique(ppi3$G1), function(x) {
        ppi2$G2[ppi2$G1 == x]
    })
    names(ppil) <- unique(ppi3$G1)

    ## Get second-order neighbors
    second_details <- Map(function(x, y) {
        second_df <- ppi2[ppi2$G1 %in% c(x, y), ]
        second_genes <- unique(c(second_df$G1, second_df$G2))
        setdiff(second_genes, c(x, y))
    }, ppil, names(ppil))

    retained <- vapply(second_details, length, integer(1)) >= min.second.neighbor.size
    ppil <- ppil[retained]
    second_details <- second_details[retained]

    message(sprintf("%d center genes detected", length(ppil)))

    ## Calculate module scores
    message("Calculating module scores...")
    score_l <- lapply(data, function(x) {
        correlation <- x$r
        v <- x$v

        scores <- do.call(rbind, Map(function(genes_in, core_gene) {
            genes_out <- second_details[[core_gene]]

            cor_in <- mean(correlation[core_gene, genes_in], na.rm = TRUE)
            cor_out <- mean(correlation[genes_in, genes_out], na.rm = TRUE)
            v_in <- mean(v[genes_in], na.rm = TRUE)

            ## Handle edge cases
            if (is.na(cor_out) || abs(cor_out) < 1e-10) cor_out <- 1e-10

            score <- if (AddModuleSize) {
                sqrt(length(genes_in)) * v_in * cor_in / cor_out
            } else {
                v_in * cor_in / cor_out
            }

            if (!is.finite(score)) score <- 0

            data.frame(
                s = score, v_in = v_in,
                r_in = cor_in, r_out = cor_out
            )
        }, ppil, names(ppil)))

        data.frame(
            Module = names(ppil),
            Rank = rank(-scores$s),
            Size = vapply(ppil, length, integer(1)),
            CI = scores$s,
            V_in = scores$v_in,
            R_in = scores$r_in,
            R_out = scores$r_out,
            row.names = NULL
        )
    })

    ## Determine number of top genes
    N <- if (percent) {
        ceiling(length(ppil) * top.p)
    } else {
        if (top.n > length(ppil)) {
            stop("top.n (", top.n, ") exceeds number of center genes (",
                 length(ppil), ")")
        }
        top.n
    }

    ## Calculate candidate statistics
    Candidate <- data.frame(
        State = names(score_l),
        CI = vapply(score_l, function(x) {
            mean(x$CI[x$Rank <= N])
        }, numeric(1)),
        row.names = NULL
    )

    ## Identify DNB genes
    Candidate_DNB_list <- lapply(score_l, function(x) {
        x$Module[x$Rank <= N]
    })
    critical_idx <- which.max(Candidate$CI)
    DNB.genes <- Candidate_DNB_list[[critical_idx]]

    ## Calculate DNB scores
    DNB.score <- data.frame(
        State = state.levels,
        DNB.score = vapply(score_l, function(x) {
            mean(x$CI[x$Module %in% DNB.genes])
        }, numeric(1))
    )

    ## Report results
    message(sprintf("Critical state: %s", Candidate$State[critical_idx]))
    message(sprintf("DNB genes identified: %d", length(DNB.genes)))
    message("Done!")

    list(
        DNB.score = DNB.score,
        DNB.genes = DNB.genes,
        CI_all = score_l,
        Gene_module = ppil,
        Candidate = Candidate,
        Cor = corl,
        V = vl,
        PPI.used = ppi2,
        first.order.genes = ppil,
        second.order.genes = second_details
    )
}

## Internal validation functions
.validate_expr <- function(expr) {
    if (!is.matrix(expr) && !is.data.frame(expr)) {
        stop("'expr' must be a matrix or data.frame", call. = FALSE)
    }
    if (is.null(rownames(expr))) {
        stop("'expr' must have row names (gene symbols)", call. = FALSE)
    }
    if (is.null(colnames(expr))) {
        stop("'expr' must have column names (sample IDs)", call. = FALSE)
    }
}

.validate_state <- function(state, sample_ids) {
    if (ncol(state) != 2) {
        stop("'state' must have exactly 2 columns: sample ID and state",
             call. = FALSE)
    }
    if (!all(state[[1]] %in% sample_ids)) {
        missing <- setdiff(state[[1]], sample_ids)
        stop("Some samples in 'state' not found in 'expr': ",
             paste(head(missing, 3), collapse = ", "),
             if (length(missing) > 3) "...", call. = FALSE)
    }
}
