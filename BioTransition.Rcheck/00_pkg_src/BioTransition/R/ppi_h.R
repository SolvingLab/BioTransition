#' Human Protein-Protein Interaction Network
#'
#' @description
#' Protein-protein interaction network for human (Homo sapiens) from STRING
#' database v12.0. This dataset contains interactions with confidence scores,
#' suitable for network-based DNB analyses including LcDNB, MDNB, LDNB, TSNMB,
#' and TSLE.
#'
#' @format A data frame with protein-protein interactions:
#' \describe{
#'   \item{G1}{Gene 1 symbol (character). First interacting gene.}
#'   \item{G2}{Gene 2 symbol (character). Second interacting gene.}
#'   \item{combined_score}{Interaction confidence score (numeric, 0-999).
#'     Higher scores indicate stronger evidence for interaction.
#'     Typical cutoffs: 400 (medium), 700 (high), 900 (highest).}
#' }
#'
#' @details
#' **Data Source**: STRING database v12.0 (https://string-db.org/)
#'
#' **Species**: Homo sapiens (NCBI Taxonomy ID: 9606)
#'
#' **Recommended Score Cutoffs**:
#' \itemize{
#'   \item 400-600: Medium confidence (broad coverage)
#'   \item 700-800: High confidence (balanced)
#'   \item 900+: Highest confidence (most reliable)
#' }
#'
#' **Usage in DNB Analysis**:
#' Use this PPI network for human single-cell or bulk RNA-seq data with
#' PPI-based DNB methods: LcDNB, MDNB, LDNB, TSNMB, TSLE.
#'
#' @source
#' STRING database v12.0
#' \url{https://string-db.org/}
#' \url{https://stringdb-downloads.org/}
#'
#' @references
#' Szklarczyk D, et al. The STRING database in 2023: protein-protein association
#' networks and functional enrichment analyses for any sequenced genome of
#' interest. Nucleic Acids Res. 2023;51(D1):D638-D646.
#'
#' @seealso
#' \code{\link{ppi_m}} for mouse PPI network
#'
#' @examples
#' # Load human PPI
#' data("ppi_h")
#'
#' # Explore the network
#' dim(ppi_h)
#' head(ppi_h)
#'
#' # Check gene coverage
#' all_genes <- unique(c(ppi_h$G1, ppi_h$G2))
#' length(all_genes)
#'
#' # Filter by confidence
#' high_conf <- ppi_h[ppi_h$combined_score >= 700, ]
#' nrow(high_conf)
#'
#' @name ppi_h
#' @docType data
#' @keywords datasets
"ppi_h"
