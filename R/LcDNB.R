#' Local Conventional DNB Analysis
#'
#' Performs local conventional dynamic network biomarker (LcDNB) analysis using
#' protein-protein interaction (PPI) networks.
#'
#' @param expr A numeric matrix with genes in rows and samples in columns.
#' @param state A data.frame with sample IDs and state labels (2 columns).
#' @param state.levels Character vector specifying the order of states.
#' @param cor.method Correlation method: "pearson", "spearman", or "kendall".
#' @param p.adjust.method P-value adjustment method. Default: "BH".
#' @param variation.method Method for variation: "sd" or "cv". Default: "sd".
#' @param min.first.neighbor.size Minimum first-order neighbors. Default: 3.
#' @param min.second.neighbor.size Minimum second-order neighbors. Default: 1.
#' @param ppi PPI network data.frame with G1, G2, combined_score columns.
#' @param min.combined.score Minimum STRING score. Default: 900.
#' @param percent Use percentage (TRUE) or absolute number (FALSE).
#' @param top.n Number of top genes when percent=FALSE. Default: 30.
#' @param top.p Proportion when percent=TRUE. Default: 0.05.
#' @param AddModuleSize Weight by module size. Default: FALSE.
#'
#' @return A list containing DNB.score, DNB.genes, CI_all, Gene_module,
#'   Candidate, Cor, V, PPI.used, first.order.genes, second.order.genes.
#'
#' @author Zaoqu Liu
#'
#' @seealso \code{\link{cDNB}}, \code{\link{LDNB}}, \code{\link{ppi_h}}
#'
#' @examples
#' \donttest{
#' # See vignette for detailed examples
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
    if (!is.matrix(expr) && !is.data.frame(expr)) {
        stop("'expr' must be a matrix or data.frame", call. = FALSE)
    }
    checkStateTwoCol(state)

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
    checkNoEmptyState(ddl)

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
            stop("top.n exceeds number of center genes")
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
