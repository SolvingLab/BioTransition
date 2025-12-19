# BioTransition

> Detect Critical Transitions in Biological Systems

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-brightgreen.svg)](https://www.r-project.org/)

## Overview

**BioTransition** identifies critical transition points and biomarker genes during biological state transitions using Dynamic Network Biomarker (DNB) theory. The package implements seven state-of-the-art DNB methods with high-performance C++ acceleration.

## Key Features

- 🔬 **7 DNB Methods**: cDNB, tDNB, LcDNB, LDNB, MDNB, TSNMB, TSLE
- ⚡ **C++ Acceleration**: Automatic optimization (2-20x speedup)
- 📊 **Multi-omics Support**: Bulk RNA-seq, scRNA-seq, spatial transcriptomics
- 🧬 **Built-in PPI Networks**: Human and mouse from STRING database
- 🎯 **Multi-group Analysis**: Support for 2+ state comparisons

## Installation

```r
# Install from GitHub
devtools::install_github("SolvingLab/BioTransition")
```

## Quick Start

```r
library(BioTransition)

# Prepare sample groups
sample_groups <- data.frame(
  sample_id = colnames(your_data),
  group = your_groups
)

# Run cDNB analysis
result <- cDNB(
  expr = your_expression_matrix,
  state = sample_groups,
  state.levels = c("Control", "Treatment")
)

# View results
result$DNB.score
result$DNB.genes
```

## Methods

| Method | Speed | PPI Required | Reference Required | Best For |
|--------|-------|--------------|-------------------|----------|
| **cDNB** | ⚡⚡⚡ | No | No | Quick screening |
| **tDNB** | ⚡⚡⚡ | No | No | Network topology |
| **LcDNB** | ⚡⚡ | Yes | No | Local networks |
| **MDNB** | ⚡ | Yes | No | Module analysis |
| **LDNB** | ⚡ | Yes | Yes | Multi-state comparison |
| **TSNMB** | ⚡⚡ | Yes | Yes | Time-series |
| **TSLE** | ⚡⚡ | Yes | Yes | Leading edge |

## References

1. Chen et al. (2012) *Sci Rep*. DOI: [10.1038/srep00342](https://doi.org/10.1038/srep00342)
2. Liu et al. (2019) *Natl Sci Rev*. DOI: [10.1093/nsr/nwy162](https://doi.org/10.1093/nsr/nwy162)
3. Li et al. (2022) *Innovation*. DOI: [10.1016/j.xinn.2022.100364](https://doi.org/10.1016/j.xinn.2022.100364)
4. Zhong et al. (2022) *J Mol Cell Biol*. DOI: [10.1093/jmcb/mjac052](https://doi.org/10.1093/jmcb/mjac052)
5. Liu et al. (2020) *Bioinformatics*. DOI: [10.1093/bioinformatics/btz758](https://doi.org/10.1093/bioinformatics/btz758)

## Author

**Zaoqu Liu**  
📧 liuzaoqu@163.com

## License

GPL (>= 3)

