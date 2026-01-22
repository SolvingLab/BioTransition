# Global variables declaration for NSE (non-standard evaluation)
# This silences R CMD check NOTEs about undefined global variables
utils::globalVariables(c(
    # dplyr/tidyverse NSE variables
    ".",
    "G1", "G2",
    "V1", "V2",
    # Data frame column names used in various functions
    "DNB.score",
    "State",
    "GI", "LI",
    "GLSE", "LSLE",
    "RLE", "RSD",
    # Data objects
    "ppi_h"
))
