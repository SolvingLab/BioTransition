#' @title Landscape DNB Analysis
#'
#' @description
#' Performs landscape dynamic network biomarker (LDNB) analysis for detecting
#' critical transitions. This method uses sample-specific perturbation networks
#' (SSPN) to identify tipping points in biological state transitions.
#'
#' @param expr A numeric matrix or data.frame with genes in rows and samples in
#'   columns. Row names should be gene symbols.
#' @param state A data.frame with exactly two columns: sample identifiers and
#'   state labels. Must include a \code{"ref"} state for reference samples.
#' @param state.levels A character vector specifying the order of states. The
#'   first level should be \code{"ref"} (reference state).
#' @param cor.method Character string specifying correlation method. One of
#'   \code{"pearson"} (default), \code{"spearman"}, or \code{"kendall"}.
#' @param p.adjust.method Character string specifying p-value adjustment method.
#'   Default: \code{"BH"}.
#' @param ppi A data.frame containing protein-protein interactions with columns
#'   \code{G1}, \code{G2}, and \code{combined_score}.
#' @param min.combined.score Numeric. Minimum STRING combined score. Default:
#'   990.
#' @param min.first.neighbor.size Integer. Minimum first-order neighbors.
#'   Default: 10.
#' @param min.second.neighbor.size Integer. Minimum second-order neighbors.
#'   Default: 1.
#' @param use.PCC.P.type Character. Type of p-value for filtering: \code{"FDR"}
#'   (default) or \code{"NP"} (nominal p-value).
#' @param use.PCC.P.cutoff Numeric. P-value cutoff. Default: 0.05.
#' @param percent Logical. Use percentage (\code{top.p}) or absolute number
#'   (\code{top.n}). Default: \code{FALSE}.
#' @param top.n Integer. Number of top genes when \code{percent = FALSE}.
#'   Default: 30.
#' @param top.p Numeric. Proportion when \code{percent = TRUE}. Default: 0.05.
#' @param nCores Integer. Number of CPU cores for parallel computation.
#'
#' @return A list containing:
#' \describe{
#'   \item{state.GI}{Data.frame with global indices per state}
#'   \item{DNB.genes}{Character vector of identified DNB genes}
#'   \item{Gene.LI}{Data.frame of genes ranked by landscape index}
#'   \item{case.GI}{Data.frame with global index per sample}
#'   \item{case.LI.list}{List of local indices per sample}
#'   \item{SSPN}{Sample-specific perturbation network results}
#' }
#'
#' @details
#' LDNB analysis requires a reference state (labeled as \code{"ref"} in the
#' state column) to construct background networks. The algorithm:
#' \enumerate{
#'   \item Constructs sample-specific perturbation networks (SSPN)
#'   \item Calculates local landscape indices for each gene
#'   \item Aggregates to global indices per sample and state
#'   \item Identifies critical state and DNB genes
#' }
#'
#' @references
#' Liu R, et al. (2019). Single-sample landscape entropy reveals the imminent
#' phase transition during disease progression. National Science Review,
#' 7(7):775-785. \doi{10.1093/nsr/nwy162}
#'
#' @author Zaoqu Liu \email{liuzaoqu@@163.com}
#'
#' @seealso \code{\link{SSPN1}}, \code{\link{LcDNB}}, \code{\link{SLE}}
#'
#' @examples
#' \donttest{
#' # LDNB requires a reference state
#' # See vignette for detailed examples
#' }
#'
#' @export
LDNB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    ppi = ppi_h,
    min.combined.score = 990,
    min.first.neighbor.size = 10,
    min.second.neighbor.size = 1,
    use.PCC.P.type = "FDR",
    use.PCC.P.cutoff = 0.05,
    percent = FALSE,
    top.n = 30,
    top.p = 0.05,
    nCores = max(1, parallel::detectCores() - 2)) {

    message("LDNB: Landscape DNB Analysis")
    message("Reference: Liu et al. (2019) Natl Sci Rev. doi:10.1093/nsr/nwy162")
    message("-------------------------------------------")

    ## Validate input
    checkStateTwoCol(state)
    colnames(state) <- c("ID", "state")
    state$state <- factor(state$state, levels = state.levels)
    state <- state[match(colnames(expr), state$ID), ]

    checkRefState(state.levels)

    ref.samples <- state$ID[state$state == "ref"]

    ## Build SSPN
    message("Building sample-specific perturbation networks...")
    SSN <- SSPN1(
        expr = expr,
        ref.samples = ref.samples,
        cor.method = cor.method,
        p.adjust.method = p.adjust.method,
        ppi = ppi,
        min.combined.score = min.combined.score,
        nCores = nCores
    )
    case.SSPN.list <- SSN$case.SSPN.list

    ## Prepare expression data
    all_genes <- unique(c(case.SSPN.list[[1]]$V1, case.SSPN.list[[1]]$V2))
    expr <- expr[rownames(expr) %in% all_genes, , drop = FALSE]
    refExpr <- expr[, ref.samples, drop = FALSE]
    refMean <- rowMeans(refExpr)

    ## Filter SSPN by p-value
    case.SSPN.list2 <- lapply(case.SSPN.list, function(x) {
        if (use.PCC.P.type == "FDR") {
            x2 <- x[x$caseFDR < use.PCC.P.cutoff, ]
        } else {
            x2 <- x[x$caseP < use.PCC.P.cutoff, ]
        }
        x3 <- x2[, c(2, 1, 3:6)]
        colnames(x3)[1:2] <- c("V1", "V2")
        rbind(x2, x3)
    })

    ## Calculate landscape indices
    message("Calculating landscape indices...")
    plan(multisession, workers = nCores)

    res.list <- furrr::future_map2(
        case.SSPN.list2, names(case.SSPN.list2),
        function(s, id) {
            .calc_landscape_index(
                s, id, expr, refMean,
                min.first.neighbor.size, min.second.neighbor.size,
                percent, top.n, top.p
            )
        },
        .progress = TRUE
    )

    plan(sequential)

    ## Extract results
    case.LI.list <- lapply(res.list, `[[`, "LI")
    case.GI <- do.call(rbind, Map(function(x, id) {
        data.frame(ID = id, GI = x$GI, stringsAsFactors = FALSE)
    }, res.list, names(res.list)))
    case.GI <- merge(state, case.GI, by = "ID")

    ## Handle failed samples
    failed_samples <- case.GI$ID[case.GI$GI == 0]
    if (length(failed_samples) > 0) {
        message(sprintf("[WARNING] %d samples failed", length(failed_samples)))
    }

    ## Calculate state-level indices
    state.GI <- aggregate(GI ~ state, data = case.GI, FUN = mean)

    ## Identify DNB genes
    critical_state <- as.character(state.GI$state[which.max(state.GI$GI)])
    valid_samples <- case.GI$ID[case.GI$state == critical_state & case.GI$GI > 0]

    if (length(valid_samples) == 0) {
        message("[WARNING] No valid samples in critical state, using all")
        valid_samples <- case.GI$ID[case.GI$GI > 0]
    }

    k <- case.LI.list[valid_samples]

    ## Find common genes and rank
    DNB.genes <- .identify_dnb_genes(k, percent, top.n, top.p)

    message(sprintf("Critical state: %s", critical_state))
    message(sprintf("DNB genes identified: %d", length(DNB.genes)))
    message("Done!")

    list(
        state.GI = state.GI,
        DNB.genes = DNB.genes,
        DNB.score = state.GI,
        Gene.LI = NULL,
        case.GI = case.GI,
        case.LI.list = case.LI.list,
        SSPN = SSN
    )
}

## Internal helper functions
.calc_landscape_index <- function(s, id, expr, refMean,
                                   min.first.neighbor.size,
                                   min.second.neighbor.size,
                                   percent, top.n, top.p) {
    if (nrow(s) == 0) {
        return(list(
            LI = data.frame(Gene = character(0), LI = numeric(0)),
            GI = 0
        ))
    }

    ## Build local networks
    sppil <- split(s$V2, s$V1)
    retained <- vapply(sppil, length, integer(1)) >= min.first.neighbor.size
    sppil <- sppil[retained]

    if (length(sppil) == 0) {
        return(list(LI = data.frame(Gene = character(0), LI = numeric(0)), GI = 0))
    }

    ## Calculate expression deviation
    sExp <- expr[, id]
    names(sExp) <- rownames(expr)
    DsExp <- abs(sExp - refMean)

    ## Calculate local indices
    LI_values <- vapply(names(sppil), function(gene) {
        neighbors <- sppil[[gene]]
        sED.in <- mean(DsExp[c(gene, neighbors)], na.rm = TRUE)
        sPCC.in <- mean(abs(s$corPert[s$V1 == gene & s$V2 %in% neighbors]),
                        na.rm = TRUE)

        s2 <- s[s$V1 %in% neighbors & !s$V2 %in% c(gene, neighbors), ]
        if (nrow(s2) >= min.second.neighbor.size) {
            sPCC.out <- mean(abs(s2$corPert), na.rm = TRUE)
            sED.in * sPCC.in / sPCC.out
        } else {
            NA_real_
        }
    }, numeric(1))

    LI_values <- LI_values[!is.na(LI_values)]

    if (length(LI_values) == 0) {
        return(list(LI = data.frame(Gene = character(0), LI = numeric(0)), GI = 0))
    }

    N <- if (percent) {
        max(1, ceiling(length(LI_values) * top.p))
    } else {
        min(top.n, length(LI_values))
    }

    GI <- mean(sort(LI_values, decreasing = TRUE)[seq_len(N)])

    list(
        LI = data.frame(Gene = names(LI_values), LI = LI_values),
        GI = GI
    )
}

.identify_dnb_genes <- function(k, percent, top.n, top.p) {
    if (length(k) == 0) return(character(0))

    all_genes <- unique(unlist(lapply(k, `[[`, "Gene")))
    if (length(all_genes) == 0) return(character(0))

    ## Calculate mean LI
    mean_LI <- vapply(all_genes, function(gene) {
        lis <- vapply(k, function(x) {
            idx <- which(x$Gene == gene)
            if (length(idx) > 0) x$LI[idx] else NA_real_
        }, numeric(1))
        mean(lis, na.rm = TRUE)
    }, numeric(1))

    mean_LI <- sort(mean_LI, decreasing = TRUE)

    N <- if (percent) {
        max(1, ceiling(length(mean_LI) * top.p))
    } else {
        min(top.n, length(mean_LI))
    }

    names(mean_LI)[seq_len(N)]
}
