#' @title Enumerate candidate gene modules by hierarchical clustering
#' @description
#' Build a dendrogram from a gene-by-gene similarity matrix (correlation or TOM)
#' and enumerate every internal-node partition as a candidate DNB module, keeping
#' only those within \code{[min.size, max.size]}.
#'
#' Module enumeration relies on \code{dendextend::partition_leaves()}, which
#' walks the tree recursively. For a very large, deeply unbalanced tree (common
#' when the input is the whole transcriptome rather than a pre-filtered gene
#' set) the recursion can exceed R's C stack and abort with a cryptic
#' \dQuote{C stack usage ... is too close to the limit} error. We intercept that
#' case and raise an actionable message instead of letting the low-level error
#' surface to the user.
#' @param simil A gene-by-gene similarity matrix (correlation or TOM). Higher
#'   values mean more similar; distance is taken as \code{1 - simil}.
#' @param min.size Minimum gene number of a module.
#' @param max.size Maximum gene number of a module. If it is greater than or
#'   equal to the gene number, it is capped at 80\% of the gene number.
#' @return A list of character vectors, one per kept module.
#' @noRd
detectModules <- function(simil, min.size, max.size) {
    n_gene <- nrow(simil)
    if (max.size >= n_gene) {
        max.size <- n_gene * 0.8
    }
    dend <- stats::as.dendrogram(stats::hclust(stats::as.dist(1 - simil)))
    modules <- tryCatch(
        dendextend::partition_leaves(dend),
        error = function(e) {
            if (grepl("C stack", conditionMessage(e), fixed = TRUE)) {
                stop(
                    "module enumeration exceeded R's C stack: the clustering ",
                    "tree for ", n_gene, " genes is too deep.\n",
                    "  DNB analysis is meant to run on a pre-filtered gene set, ",
                    "not the whole transcriptome.\n",
                    "  Reduce the input first, e.g. keep the top 1000-3000 most ",
                    "variable or differentially\n",
                    "  expressed genes, then rerun. (Current gene number: ",
                    n_gene, ")",
                    call. = FALSE
                )
            }
            stop(e)
        }
    )
    member_nums <- vapply(modules, length, integer(1))
    modules[member_nums >= min.size & member_nums <= max.size]
}

#' @title Warn when the input gene number is large
#' @description
#' DNB methods are designed for a pre-filtered gene set. A one-line heads-up
#' when the input is large nudges the user to filter before a slow run or a C
#' stack overflow during module detection.
#' @param n_gene Number of input genes (\code{nrow(expr)}).
#' @param threshold Gene number above which to emit the note.
#' @return Invisibly \code{NULL}; called for its message side effect.
#' @noRd
warnLargeGeneNumber <- function(n_gene, threshold = 5000) {
    if (n_gene > threshold) {
        message("+++ Note: input has ", n_gene, " genes. DNB analysis is ",
                "usually run on a pre-filtered set")
        message("+++ (e.g. top 1000-3000 variable or DE genes); very large ",
                "gene numbers slow clustering")
        message("+++ and may exceed R's C stack during module detection.")
        message("")
    }
    invisible(NULL)
}

#' @title Score candidate gene modules with the DNB composite index
#' @description
#' Compute the DNB composite index (CI) for each candidate module under the
#' \code{|PCC|} convention of the original DNB theory (Chen et al., 2012):
#' \deqn{CI = mean(SD_in) * mean(|PCC_in|) / mean(|PCC_out|)}
#'
#' Performance note: the score is dominated by \code{mean(|PCC_out|)}, the mean
#' absolute similarity between a module and every gene outside it. Computing it
#' naively per module is \eqn{O(n_{in} \cdot N)} and forces a fresh
#' module-vs-outside submatrix copy each time. We avoid that with a complement
#' identity: with \code{A = |simil|} and its per-gene row sums precomputed once,
#' \deqn{sum|PCC(in, out)| = sum(rowSums(A)[in]) - sum(A[in, in])}
#' so each module only touches its own \eqn{n_{in} \times n_{in}} block. This
#' is far faster than both the per-module R loop and the bundled
#' \code{fast_module_score_cpp} (whose full-matrix copy and linear is-outside
#' scan make it ~50x slower than this vectorised form in practice).
#' @param simil A gene-by-gene similarity matrix (correlation or TOM) carrying
#'   gene names on both rows and columns.
#' @param variation A named numeric vector of per-gene variation.
#' @param modules A list of character vectors, each one candidate module.
#' @param add.size Whether to weight CI by \code{sqrt(module size)}.
#' @return A data.frame with columns Module, Size, CI, V_in, R_in, R_out.
#' @noRd
scoreModules <- function(simil, variation, modules, add.size = FALSE) {
    genes <- rownames(simil)
    variation <- variation[genes]
    n_total <- length(genes)
    abs_simil <- abs(simil)
    row_sum <- rowSums(abs_simil)
    gidx <- stats::setNames(seq_along(genes), genes)
    sizes <- vapply(modules, length, integer(1))

    stat <- vapply(modules, function(mod) {
        mi <- gidx[mod]
        n_in <- length(mi)
        sub <- abs_simil[mi, mi, drop = FALSE]
        sum_in_all <- sum(sub)
        # within-module mean |PCC|, excluding the self-correlation diagonal
        sum_in_off <- sum_in_all - sum(diag(sub))
        pairs_in <- n_in * (n_in - 1)
        r_in <- if (pairs_in > 0) sum_in_off / pairs_in else 0
        # module-vs-outside mean |PCC| via the complement identity
        sum_in_out <- sum(row_sum[mi]) - sum_in_all
        pairs_out <- n_in * (n_total - n_in)
        r_out <- if (pairs_out > 0) sum_in_out / pairs_out else 0
        v_in <- mean(variation[mi])
        ci <- if (add.size) sqrt(n_in) * v_in * r_in / r_out else v_in * r_in / r_out
        c(ci, v_in, r_in, r_out)
    }, numeric(4))

    data.frame(
        Module = seq_along(modules),
        Size = sizes,
        CI = stat[1, ],
        V_in = stat[2, ],
        R_in = stat[3, ],
        R_out = stat[4, ]
    )
}
