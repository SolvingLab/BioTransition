#' @importFrom Rcpp sourceCpp
#' @useDynLib BioTransition, .registration = TRUE
NULL

# Global variables declaration for NSE (non-standard evaluation)
# This silences R CMD check NOTEs about undefined global variables
utils::globalVariables(c(
    # dplyr/tidyverse NSE variables
    ".",
    "G1", "G2",
    "V1", "V2",
    # Data frame column names
    "DNB.score", "State", "state",
    "GI", "LI", "CI",
    "GLSE", "LSLE",
    "RLE", "RSD",
    "Rank", "Module",
    "caseFDR", "caseP", "corPert",
    # Data objects
    "ppi_h", "ppi_m"
))

#' @keywords internal
.onLoad <- function(libname, pkgname) {
    # Set default options if needed
    invisible(NULL)
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "BioTransition v", utils::packageVersion(pkgname), "\n",
        "Dynamic Network Biomarker Analysis for Critical Transition Detection\n",
        "GitHub: https://github.com/SolvingLab/BioTransition"
    )
}
