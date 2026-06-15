#' @title Time-series network module biomarker analysis
#' @description
#' Performing time-series network module biomarker analysis.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param state A time-series dataframe with two columns, the first is the sample names and the second is the group or time point information.
#' @param state.levels A vector for state sequence.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param min.first.neighbor.size Minimum size of first order genes of a specific center gene.
#' @param ppi Protein-protein interaction network; background network.
#' @param min.combined.score Minimum combined score for determining protein-protein interaction.
#' @param percent Whether to use Percent to determine the number of DNB genes.
#' @param top.n Only percent = FALSE takes effect. Center genes with top (number) DNB score were defined as DNB genes.
#' @param top.p Only percent = TRUE takes effect. Center genes with top (percent) DNB score were defined as DNB genes.
#' @import magrittr
#' @export
TSNMB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    ppi = ppi_h,
    min.combined.score = 900,
    min.first.neighbor.size = 3,
    percent = TRUE,
    top.n = 30,
    top.p = 0.05) {
  message("+++ ThiS method was modified by Zaoqu Liu on the basis of sNMB!")
  message("")

  message("+++ Performing time-series network module biomarker analysis...")
  message("")
  checkStateTwoCol(state)
  colnames(state) <- c("ID", "state")
  state$state <- factor(state$state, state.levels)
  state <- state[match(colnames(expr), state$ID), ]
  checkRefState(state.levels)
  ddl <- purrr::map(state.levels, \(x){
    expr[, state$ID[state$state == x]]
  })
  names(ddl) <- state.levels

  ppi2 <- ppi[ppi$G1 %in% rownames(expr) & ppi$G2 %in% rownames(expr) & ppi$combined_score >= min.combined.score, ]

  expr <- expr[unique(c(ppi2$G1)), ]
  message("+++ ", nrow(expr), " intersected genes between expression matrix and PPI network...")
  message("")

  ppi3 <- ppi2[ppi2$G1 %in% names(table(ppi2$G1))[table(ppi2$G1) >= min.first.neighbor.size], ]

  ppil <- purrr::map(unique(ppi3$G1), \(x){
    g <- ppi2[ppi2$G1 == x, ]
    return(g$G2)
  })
  names(ppil) <- unique(ppi3$G1)

  message("+++ ", length(ppil), " centre genes detected in this dataset...")
  message("")
  message("+++ Calculating background network using reference samples...")
  message("")

  refExpr <- expr[, state$ID[state$state == "ref"]]

  refSD <- apply(refExpr, 1, stats::sd)

  refLNSD <- refSD[names(ppil)]

  refCor <- suppressWarnings(CorandPval(t(refExpr),
    method = cor.method,
    p.adjust.method = p.adjust.method
  ))$r

  refNet <- purrr::map2(names(ppil), ppil, ~ data.frame(Gene = .y, PCC = refCor[.y, .x], SD = refSD[.y], row.names = NULL))
  names(refNet) <- names(ppil)

  state_case <- state[state$state != "ref", ]

  if (percent == FALSE) {
    if (top.n > length(ppil)) {
      stop(paste0(
        "\nThe number of candidata core genes is ",
        length(ppil), ", which is larger than the top.n you inputted!\n",
        "Please change another number!"
      ))
    }
    N <- top.n
  } else {
    N <- ceiling(length(ppil) * top.p)
  }

  tmp_fun <- function(G) {
    caseExpr <- cbind(refExpr, case = expr[, state_case$ID[state_case$state == G]])
    caseSD <- apply(caseExpr, 1, stats::sd)
    caseLNSD <- caseSD[names(ppil)]
    DLNSD <- abs(caseLNSD - refLNSD)
    caseCor <- suppressWarnings(CorandPval(t(caseExpr), method = cor.method, p.adjust.method = p.adjust.method))$r
    caseNet <- purrr::map2(names(ppil), ppil, ~ data.frame(Gene = .y, PCC = caseCor[.y, .x], SD = caseSD[.y], row.names = NULL))
    names(caseNet) <- names(ppil)
    DNet <- purrr::map2(refNet, caseNet, \(x, y){
      data.frame(Gene = x$Gene, DPCC = abs(x$PCC - y$PCC), DSD = abs(x$SD - y$SD))
    })
    LI <- purrr::map2_vec(DNet, DLNSD, \(x, y){
      ((sum(x$DSD) + y) / (nrow(x) + 1)) * (sum(x$DPCC) / nrow(x))
    })
    LI <- data.frame(Gene = names(DLNSD), LI = LI, row.names = NULL) %>%
      dplyr::arrange(dplyr::desc(LI))
    GI <- sum(LI$LI[seq_len(N)]) / N
    return(list(DNet = DNet, LI = LI, GI = GI))
  }

  message("+++ Calculating local sNMB score for each state...")
  message("")
  state.res <- purrr::map(state.levels[-1], ~ tmp_fun(.x), .progress = TRUE)
  names(state.res) <- state.levels[-1]
  state.LI.list <- purrr::map(state.res, ~ .x$LI)

  message("")
  message("+++ Calculating global sNMB score for each state...")
  message("")
  state.GI <- data.frame(state = names(state.res), GI = purrr::map_vec(state.res, ~ .x$GI), row.names = NULL)

  state.DNet <- purrr::map(state.res, ~ .x$DNet)

  signaling.gene <- state.LI.list[[as.character(state.GI$state[which.max(state.GI$GI)])]] %>%
    dplyr::arrange(dplyr::desc(LI)) %>%
    dplyr::top_n(n = N, wt = LI) %>%
    .[, 1]

  message("+++ ", state.GI$state[which.max(state.GI$GI)], " is the critical state...")
  message("+++ ", length(signaling.gene), " signaling genes involved in the critical state...")
  message("")
  message("+++ Done!")

  return(list(
    state.GI = state.GI,
    state.LI.list = state.LI.list,
    signaling.gene = signaling.gene,
    state.DNet = state.DNet,
    refNet = refNet,
    refLNSD = refLNSD,
    Gene_module = ppil,
    PPI.used = ppi2,
    DNB.genes = signaling.gene,
    DNB.score = state.GI
  ))
}
