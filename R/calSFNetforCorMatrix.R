#' @title Calculating scale-free network based on correlation matrix
#' @description
#'  Calculate the adjacency matrix, topological matrix, and identify the optimal soft threshold for defining scale-free network.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param cordata a correlation matrix.
#' @param power.vec a vector of soft thresholding powers for which the scale free topology fit indices are to be calculated.
#' @param network.type network type. Allowed values are (unique abbreviations of) "unsigned" or "signed".
#' @param RsquaredCut desired minimum scale free topology fitting index R^2.
#' @noRd
calSFNetforCorMatrix <- function(
    cordata,
    power.vec = c(c(1:10), seq(from = 12, to = 20, by = 2)),
    network.type = "signed",
    RsquaredCut = 0.85) {
  tmpfun <- function(cordata, power, network.type = "signed") {
    if (network.type == "signed") {
      ADJ <- abs(0.5 + 0.5 * cordata)^power
    }
    if (network.type == "unsigned") {
      ADJ <- abs(cordata)^power
    }
    k <- apply(ADJ, 2, sum) - 1 # Connectivity
    cut1 <- cut(k, 10)
    binned.k <- tapply(k, cut1, mean) # Mean Connectivity for each Bin
    freq1 <- tapply(k, cut1, length) / length(k)
    xx <- as.vector(log10(binned.k))
    lm1 <- stats::lm(as.numeric(log10(freq1 + .000000001)) ~ xx)
    return(data.frame(
      Power = power,
      SFT.R.sq = summary(lm1)$r.squared,
      slope = lm1$coefficients[[2]],
      mean.k = mean(k)
    ))
  }
  pickSotfThres <- purrr::map_df(power.vec, \(x)tmpfun(cordata,
    power = x,
    network.type = network.type
  ))
  pick_power <- ifelse(max(pickSotfThres$SFT.R.sq) >= RsquaredCut,
    which(pickSotfThres$SFT.R.sq >= RsquaredCut)[1],
    pickSotfThres$Power[which.max(pickSotfThres$SFT.R.sq)]
  )

  if (network.type == "signed") {
    ADJ <- abs(0.5 + 0.5 * cordata)^pick_power
  }
  if (network.type == "unsigned") {
    ADJ <- abs(cordata)^pick_power
  }
  diag(ADJ) <- 0
  ADJ[is.na(ADJ)] <- 0

  # Topological overlap (Zhang & Horvath 2005), computed on the full matrix so
  # the denominator can be guarded and the result clamped to a valid [0, 1]
  # similarity:  TOM_ij = (L_ij + a_ij) / (min(k_i, k_j) + 1 - a_ij), L = ADJ^2
  kk <- colSums(ADJ) # Connectivity for each gene
  numTOM <- ADJ %*% ADJ + ADJ
  denomTOM <- outer(kk, kk, pmin) + 1 - ADJ
  TOMmatrix <- numTOM / denomTOM
  # Guard the degenerate denominator (isolated genes / a_ij -> 1), then keep a
  # proper similarity with unit self-overlap.
  TOMmatrix[!is.finite(TOMmatrix)] <- 0
  TOMmatrix <- pmin(pmax(TOMmatrix, 0), 1)
  diag(TOMmatrix) <- 1

  return(list(
    pickSotfThres = pickSotfThres,
    pick_power = pick_power,
    Connectivity = kk,
    ADJmatrix = ADJ,
    TOMmatrix = TOMmatrix
  ))
}
