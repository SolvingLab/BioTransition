#' Mouse Protein-Protein Interaction Network
#'
#' @description
#' Protein-protein interaction network for mouse (Mus musculus) from STRING database v12.0.
#' This dataset contains interactions with confidence scores, suitable for network-based
#' DNB analyses including LcDNB, MDNB, LDNB, TSNMB, and TSLE.
#'
#' @format A data frame with 12,684,354 protein-protein interactions:
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
#' **Species**: Mus musculus (NCBI Taxonomy ID: 10090)
#'
#' **Statistics**:
#' \itemize{
#'   \item Total interactions: 12,684,354
#'   \item Unique genes: 21,645
#'   \item Score range: 150-999
#' }
#'
#' **Recommended Score Cutoffs**:
#' \itemize{
#'   \item 400-600: Medium confidence (broad coverage)
#'   \item 700-800: High confidence (balanced)
#'   \item 900+: Highest confidence (most reliable)
#' }
#'
#' **Usage in DNB Analysis**:
#' Use this PPI network for mouse single-cell or bulk RNA-seq data with
#' PPI-based DNB methods: LcDNB, MDNB, LDNB, TSNMB, TSLE.
#'
#' @source
#' STRING database v12.0
#' \url{https://string-db.org/}
#' \url{https://stringdb-downloads.org/}
#'
#' @references
#' Szklarczyk D, et al. The STRING database in 2023: protein-protein association
#' networks and functional enrichment analyses for any sequenced genome of interest.
#' Nucleic Acids Res. 2023;51(D1):D638-D646.
#'
#' @seealso
#' \code{\link{ppi_h}} for human PPI network
#'
#' @examples
#' # Load mouse PPI
#' data("ppi_m")
#' 
#' # Explore the network
#' dim(ppi_m)
#' head(ppi_m)
#' 
#' # Check gene coverage
#' all_genes <- unique(c(ppi_m$G1, ppi_m$G2))
#' length(all_genes)  # 21,645 genes
#' 
#' # Filter by confidence
#' high_conf <- ppi_m[ppi_m$combined_score >= 700, ]
#' nrow(high_conf)
#' 
#' \dontrun{
#' # Use in DNB analysis
#' result <- MDNB(
#'   expr = mouse_expr,
#'   state = sample_groups,
#'   state.levels = c("Control", "Treatment"),
#'   ppi = ppi_m,
#'   min.combined.score = 700
#' )
#' }
#'
#' @name ppi_m
#' @docType data
#' @keywords datasets
"ppi_m"

