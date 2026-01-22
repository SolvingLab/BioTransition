#' @title BioTransition: Dynamic Network Biomarker Analysis for Critical
#'   Transition Detection
#'
#' @description
#' A comprehensive toolkit for detecting critical transitions and identifying
#' dynamic network biomarkers (DNB) in biological systems. Critical transitions
#' are abrupt shifts between alternative stable states that are ubiquitous in
#' complex biological processes including disease progression, cellular
#' differentiation, and developmental transitions.
#'
#' @details
#' The package implements seven complementary DNB methodologies:
#'
#' \strong{PPI-independent methods:}
#' \itemize{
#'   \item \code{\link{cDNB}}: Conventional DNB based on the original DNB theory
#'   \item \code{\link{tDNB}}: Topological DNB using scale-free network topology
#'     (novel method developed in this package)
#' }
#'
#' \strong{PPI-dependent methods:}
#' \itemize{
#'   \item \code{\link{LcDNB}}: Local DNB leveraging PPI networks
#'   \item \code{\link{LDNB}}: Landscape DNB for state transition quantification
#'   \item \code{\link{MDNB}}: Module-based DNB for modular network analysis
#'   \item \code{\link{TSNMB}}: Time-series network module biomarker
#'   \item \code{\link{TSLE}}: Time-series leading edge analysis
#' }
#'
#' \strong{Key features:}
#' \itemize{
#'   \item C++ acceleration via Rcpp for computational efficiency
#'   \item Support for bulk RNA-seq, scRNA-seq, and spatial transcriptomics
#'   \item Built-in PPI networks for human (\code{\link{ppi_h}}) and mouse
#'     (\code{\link{ppi_m}})
#'   \item Unified API across all methods
#' }
#'
#' @section Theoretical Background:
#' The DNB theory, proposed by Chen et al. (2012), provides a model-free
#' approach to detect early warning signals of critical transitions. The core
#' principle is that a dominant group of molecules (DNB module) exhibits three
#' key characteristics near the tipping point:
#' \enumerate{
#'   \item Increased variation: SD of DNB members increases dramatically
#'   \item Enhanced correlation: PCC among DNB members strengthens
#'   \item Decreased correlation: PCC between DNB and non-DNB members weakens
#' }
#'
#' The composite index (CI) quantifies the transition signal:
#' \deqn{CI = \frac{\overline{SD_{in}} \cdot \overline{|PCC_{in}|}}
#'   {\overline{|PCC_{out}|}}}
#'
#' @references
#' Chen L, Liu R, Liu ZP, Li M, Aihara K. (2012). Detecting early-warning
#' signals for sudden deterioration of complex diseases by dynamical network
#' biomarkers. Scientific Reports, 2:342. \doi{10.1038/srep00342}
#'
#' Liu R, Chen P, Chen L. (2019). Single-sample landscape entropy reveals the
#' imminent phase transition during disease progression. National Science
#' Review, 7(7):1175-1185. \doi{10.1093/nsr/nwy162}
#'
#' Li M, Zeng T, Liu R, Chen L. (2022). Detecting tissue-specific early warning
#' signals for complex diseases based on dynamical network biomarkers. The
#' Innovation, 3(5):100364. \doi{10.1016/j.xinn.2022.100364}
#'
#' Zhong J, Han C, Wang Y, Chen P, Liu R. (2022). Identifying the critical
#' state of complex diseases by the novel concept of sample-specific network
#' module biomarker. Journal of Molecular Cell Biology, 14(6):mjac052.
#' \doi{10.1093/jmcb/mjac052}
#'
#' Liu X, Chang X, Liu R, Yu X, Chen L, Aihara K. (2020). Quantifying critical
#' states of complex diseases using single-sample dynamic network biomarkers.
#' Bioinformatics, 36(4):1068-1074. \doi{10.1093/bioinformatics/btz758}
#'
#' @author Zaoqu Liu \email{liuzaoqu@@163.com}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/SolvingLab/BioTransition}
#'   \item \url{https://zaoqu-liu.r-universe.dev/BioTransition}
#'   \item Report bugs at
#'     \url{https://github.com/SolvingLab/BioTransition/issues}
#' }
#'
#' @docType package
#' @name BioTransition-package
#' @aliases BioTransition
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats aggregate cor p.adjust sd setNames var
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @importFrom methods is
#' @import future
#' @import magrittr
#' @useDynLib BioTransition, .registration = TRUE
## usethis namespace: end
NULL
