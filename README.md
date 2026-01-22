# BioTransition <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/BioTransition)](https://zaoqu-liu.r-universe.dev/BioTransition)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-brightgreen.svg)](https://www.r-project.org/)
<!-- badges: end -->

## Overview

**BioTransition** is a comprehensive R package for detecting critical transitions in biological systems using Dynamic Network Biomarker (DNB) theory. Critical transitions—abrupt shifts between alternative stable states—are ubiquitous in complex biological processes, including disease onset, cellular differentiation, and developmental transitions. Identifying the pre-transition state (tipping point) before an irreversible shift occurs is crucial for early warning and therapeutic intervention.

This package implements seven complementary DNB methodologies, including **tDNB** (topological DNB), a novel algorithm that leverages network topology and scale-free properties for robust critical transition detection.

## Theoretical Background

The DNB theory, proposed by Chen et al. (2012), provides a model-free approach to detect early warning signals of critical transitions. The core principle is that a dominant group of molecules (DNB module) exhibits three key characteristics near the tipping point:

1. **Increased variation**: Standard deviation (SD) of DNB members increases dramatically
2. **Enhanced correlation**: Pearson correlation coefficient (PCC) among DNB members strengthens
3. **Decreased correlation**: PCC between DNB members and non-DNB molecules weakens

The composite index (CI) quantifies the transition signal:

$$CI = \frac{\overline{SD_{in}} \cdot \overline{|PCC_{in}|}}{\overline{|PCC_{out}|}}$$

## Implemented Methods

| Method | Description | Network | Reference |
|:-------|:------------|:-------:|:----------|
| **cDNB** | Conventional DNB based on correlation screening | — | Chen et al., 2012 |
| **tDNB** | Topological DNB using scale-free network topology | — | Liu Z. *(this package)* |
| **LcDNB** | Local DNB with protein-protein interaction constraints | PPI | — |
| **LDNB** | Landscape DNB for state transition quantification | PPI | Liu et al., 2019 |
| **MDNB** | Module-based DNB for modular network analysis | PPI | Li et al., 2022 |
| **TSNMB** | Time-series network module biomarker | PPI | Zhong et al., 2022 |
| **TSLE** | Time-series leading edge analysis | PPI | Liu et al., 2020 |

## Installation

```r
# Install from R-universe (recommended)
install.packages("BioTransition", repos = "https://zaoqu-liu.r-universe.dev")

# Or install from GitHub
# install.packages("devtools")
devtools::install_github("SolvingLab/BioTransition")
```

## Quick Start

```r
library(BioTransition)

# Prepare sample annotation
sample_info <- data.frame(
  sample_id = colnames(expression_matrix),
  state = state_labels
)

# Run tDNB analysis (topological DNB)
result <- tDNB(
  expr = expression_matrix,
  state = sample_info,
  state.levels = c("Normal", "Pre-disease", "Disease")
)

# Examine results
result$DNB.score    # Composite index across states
result$DNB.genes    # Identified DNB module genes
result$Candidate    # Candidate modules per state
```

## Data Compatibility

BioTransition supports multiple omics data types:

- **Bulk RNA-seq**: Population-level transcriptomics
- **Single-cell RNA-seq**: Cell-type specific analysis  
- **Spatial transcriptomics**: Spatially resolved expression data

## Built-in Resources

The package includes curated protein-protein interaction (PPI) networks:

- `ppi_h`: Human PPI network from STRING database
- `ppi_m`: Mouse PPI network from STRING database

## Performance

Core computational routines are implemented in C++ via Rcpp, providing 2–20× speedup compared to pure R implementations for correlation analysis and module scoring.

## References

1. Chen L, Liu R, Liu ZP, Li M, Aihara K. (2012). Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. *Scientific Reports*, 2:342. doi: [10.1038/srep00342](https://doi.org/10.1038/srep00342)

2. Liu R, Chen P, Chen L. (2019). Single-sample landscape entropy reveals the imminent phase transition during disease progression. *National Science Review*, 7(7):1175-1185. doi: [10.1093/nsr/nwy162](https://doi.org/10.1093/nsr/nwy162)

3. Li M, Zeng T, Liu R, Chen L. (2022). Detecting tissue-specific early warning signals for complex diseases based on dynamical network biomarkers: study of type 2 diabetes by cross-tissue analysis. *The Innovation*, 3(5):100364. doi: [10.1016/j.xinn.2022.100364](https://doi.org/10.1016/j.xinn.2022.100364)

4. Zhong J, Han C, Wang Y, Chen P, Liu R. (2022). Identifying the critical state of complex diseases by the novel concept of sample-specific network module biomarker. *Journal of Molecular Cell Biology*, 14(6):mjac052. doi: [10.1093/jmcb/mjac052](https://doi.org/10.1093/jmcb/mjac052)

5. Liu X, Chang X, Liu R, Yu X, Chen L, Aihara K. (2020). Quantifying critical states of complex diseases using single-sample dynamic network biomarkers. *Bioinformatics*, 36(4):1068-1074. doi: [10.1093/bioinformatics/btz758](https://doi.org/10.1093/bioinformatics/btz758)

## Citation

If you use BioTransition in your research, please cite:

```
Liu Z (2026). BioTransition: Dynamic Network Biomarker Analysis for Critical 
Transition Detection. R package version 2.0.0. 
https://github.com/SolvingLab/BioTransition
```

## Author

**Zaoqu Liu**, MD, PhD  
Chinese Academy of Medical Sciences and Peking Union Medical College  
Email: liuzaoqu@163.com  
ORCID: [0000-0002-0452-742X](https://orcid.org/0000-0002-0452-742X)

## License

GPL (>= 3)
