# BioTransition 2.1.0

## Major changes

* **Composite index unified to the |PCC| convention.** `cDNB()` and `tDNB()`
  now score modules with mean *absolute* correlation, matching the original DNB
  definition of Chen et al. (2012): `CI = mean(SD_in) * mean(|PCC_in|) /
  mean(|PCC_out|)`. The previous signed implementation could return negative,
  unstable CI values. **Numerical results from `cDNB()`/`tDNB()` will differ
  from 2.0.0** (and are more stable); `tDNB()` is unaffected in practice because
  its TOM is non-negative.

## Performance

* Module scoring in `cDNB()`/`tDNB()` rewritten with a complement identity
  (`sum|PCC(in,out)| = sum(rowSums|A|[in]) - sum|A[in,in]|`) and one-time
  precomputation of `|A|` and its row sums, giving ~19x speedup on realistic
  loads (3000 genes x 6000 modules: 34s -> 1.8s).

## Robustness & user experience

* Large gene sets that overflow the C stack during hierarchical clustering now
  raise an actionable message (pre-filter to the top variable/DE genes) instead
  of a cryptic "C stack usage ... too close to the limit".
* Unified, actionable input validation across all methods: empty states,
  a missing `"ref"` state, and empty PPI/case sets are reported clearly.
* `nCores` default no longer goes negative on machines with few cores
  (`SSPN1`/`SSPN2`/`SLE`/`sNMB`): now `max(1, detectCores() - 1)`.
* Progress messages unified to `message()` so they can be silenced with
  `suppressMessages()`.
* All methods now expose unified `DNB.genes` and `DNB.score` fields (existing
  method-specific field names are kept as aliases).

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
