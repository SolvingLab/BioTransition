#' @title Conventional DNB analysis
#' @description
#' Performing conventional dynamic network biomarker analysis based on the
#' original DNB theory proposed by Chen et al. (2012).
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param state A time-series dataframe with two columns, the first is the
#'   sample names and the second is the group or time point information.
#' @param state.levels A vector for state sequence.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be
#'   abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'   "fdr", "none")
#' @param variation.method specifies the method for calculating gene variation.
#'   sd or cv.
#' @param min.size Minimum gene number of gene modules.
#' @param max.size Maximum gene number of gene modules.
#' @param AddModuleSize Whether to consider gene module size when calculating
#'   DNB score.
#' @return A list containing:
#' \describe{
#'   \item{DNB.score}{A data.frame with composite index (CI) values for each
#'     state, including V_in (mean variation), R_in (mean correlation within
#'     DNB), and R_out (mean correlation with non-DNB genes)}
#'   \item{DNB.genes}{Character vector of genes in the identified DNB module}
#'   \item{CI_all}{List of data.frames with CI scores for all candidate modules
#'     in each state}
#'   \item{Gene_module}{List of gene modules detected in each state}
#'   \item{Candidate}{Data.frame summarizing the best candidate module per
#'     state}
#'   \item{Cor}{List of correlation matrices for each state}
#'   \item{V}{List of gene variation values for each state}
#' }
#' @references
#' Chen L, Liu R, Liu ZP, Li M, Aihara K. (2012). Detecting early-warning
#' signals for sudden deterioration of complex diseases by dynamical network
#' biomarkers. Scientific Reports, 2:342. \doi{10.1038/srep00342}
#'
#' @seealso
#' \code{\link{tDNB}} for topological DNB,
#' \code{\link{LcDNB}} for local DNB with PPI,
#' \code{\link{MDNB}} for module-based DNB
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' n_genes <- 100
#' n_samples <- 15
#'
#' expr <- matrix(
#'   rnorm(n_genes * n_samples, mean = 10, sd = 2),
#'   nrow = n_genes, ncol = n_samples
#' )
#' rownames(expr) <- paste0("Gene", seq_len(n_genes))
#' colnames(expr) <- paste0("Sample", seq_len(n_samples))
#'
#' state <- data.frame(
#'   sample_id = colnames(expr),
#'   state = rep(c("A", "B", "C"), each = 5)
#' )
#'
#' \donttest{
#' result <- cDNB(
#'   expr = expr,
#'   state = state,
#'   state.levels = c("A", "B", "C"),
#'   min.size = 5,
#'   max.size = 50
#' )
#' result$DNB.score
#' }
#' @export
cDNB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    variation.method = "sd",
    min.size = 10,
    max.size = 2000,
    AddModuleSize = FALSE) {
    message("Reference: Chen L et al. Detecting early-warning signals for ",
            "sudden deterioration of complex diseases by dynamical network ",
            "biomarkers. Sci Rep. 2012;2:342. doi: 10.1038/srep00342")
    message("")

    message("+++ Performing conventional dynamic network biomarker analysis...")
    message("")
    checkStateTwoCol(state)
    colnames(state) <- c("ID", "state")
    state$state <- factor(state$state, state.levels)
    state <- state[match(colnames(expr), state$ID), ]
    ddl <- purrr::map(state.levels, \(x) {
        expr[, state$ID[state$state == x]]
    })
    names(ddl) <- state.levels
    checkNoEmptyState(ddl)

    message("+++ Calculating gene correlations for each group or time point...")
    message("")
    corl <- suppressWarnings(purrr::map(ddl, ~ CorandPval(t(.x),
        method = cor.method,
        p.adjust.method = p.adjust.method
    )))

    message("+++ Calculating gene variations for each group or time point...")
    message("")
    if (variation.method == "sd") {
        vl <- purrr::map(ddl, ~ apply(.x, 1, sd))
    }
    if (variation.method == "cv") {
        vl <- purrr::map(ddl, ~ apply(.x, 1, \(a) {
            sd(a) / mean(a)
        }))
    }

    data <- purrr::map2(corl, vl, \(x, y) {
        list(r = x$r, v = y)
    })


    message("+++ Performing Hierarchical clustering for all genes ...")
    message("")
    warnLargeGeneNumber(nrow(expr))

    modulel <- purrr::map(corl, ~ detectModules(.x$r, min.size, max.size))

    message("+++ Number of gene modules detected in different groups or times:")
    module_counts <- vapply(modulel, length, integer(1))
    message(paste0(paste0("+++ ", names(modulel), " = ", module_counts,
                          collapse = "\n")))
    message("")
    message("+++ Calculating module scores for each group or time point...")
    message("")

    score_l <- purrr::map2(data, modulel,
        ~ scoreModules(.x$r, .x$v, .y, AddModuleSize))

    Candidate <- data.frame(
        State = names(score_l),
        Module = purrr::map_vec(score_l, ~ which.max(.x$CI)),
        CI = purrr::map_vec(score_l, ~ max(.x$CI)),
        V_in = purrr::map_vec(score_l, ~ .x$V_in[which.max(.x$CI)]),
        R_in = purrr::map_vec(score_l, ~ .x$R_in[which.max(.x$CI)]),
        R_out = purrr::map_vec(score_l, ~ .x$R_out[which.max(.x$CI)]),
        row.names = NULL
    )

    critical_state <- Candidate$State[which.max(Candidate$CI)]
    critical_module <- Candidate$Module[which.max(Candidate$CI)]
    DNB.genes <- modulel[[critical_state]][[critical_module]]

    DNB.score <- purrr::map_df(data, \(x) {
        s <- scoreModules(x$r, x$v, list(DNB.genes), AddModuleSize)
        data.frame(CI = s$CI, V_in = s$V_in, R_in = s$R_in, R_out = s$R_out)
    })
    DNB.score <- cbind(State = state.levels, DNB.score)

    message(paste0("+++ ", critical_state, " is the critical state..."))
    message(paste0("+++ ", length(DNB.genes),
                   " genes (DNB module) involved in the critical state..."))
    message("")
    message("+++ Done!")

    return(list(
        DNB.score = DNB.score, DNB.genes = DNB.genes,
        CI_all = score_l, Gene_module = modulel,
        Candidate = Candidate, Cor = corl, V = vl
    ))
}
