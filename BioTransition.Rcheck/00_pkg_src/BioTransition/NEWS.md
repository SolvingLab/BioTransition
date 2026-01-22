# BioTransition 2.0.0

## New Features

* Added `tDNB()` function implementing topological DNB analysis using 
  scale-free network topology
* Implemented C++ acceleration for correlation analysis and module scoring 
  via Rcpp (2-20x speedup)
* Added built-in protein-protein interaction networks for human (`ppi_h`) 
  and mouse (`ppi_m`) from STRING database
* New helper functions: `fast_cor_cpp()`, `fast_cor_pval_cpp()`, 
  `fast_module_score_cpp()`, `fast_bh_adjust()`

## Improvements

* Unified API across all seven DNB methods
* Enhanced documentation with comprehensive vignette
* Added support for single-cell RNA-seq and spatial transcriptomics data
* Improved memory efficiency for large-scale datasets

## Bug Fixes
 
* Fixed correlation calculation for samples with zero variance
* Resolved edge cases in module detection with small sample sizes

# BioTransition 1.0.0

## Initial Release

* Implemented six DNB methods: cDNB, LcDNB, LDNB, MDNB, TSNMB, TSLE
* Support for bulk RNA-seq data analysis
* Basic documentation and examples
