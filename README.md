# EGFRmSig

LUAD EGFR mSig Score Calculation


## Description

_EGFRmSig_ offers a centroid-based method to predict an EGFR mutation signature (mSig) from gene expression data in lung adenocarcinoma (LUAD). The package implements predefined EGFR wild-type (WT) and mutant (mt) centroids and applies a nearest-centroid approach to classify samples as EGFR mSig WT-like or EGFR mSig mutant-like.

See the package vignette for an application workflow: [_Tutorial_](https://miyeonyeon.github.io/bioc-vignettes/EGFRmSig_intro.html).


## Installation

```r
if (!requireNamespace("devtools", quietly=TRUE)) {
  install.packages("devtools")
}
devtools::install_github("minjeongjaykim/EGFRmSig", dependencies=TRUE)
```


## Reference
Kim M, Lamlertthon W, Jo H, et al. RAS signaling in lung adenocarcinoma is defined by lineage context and DUSP4 loss, _JCI insight_ 2026. (In Press)
Kim M, Lamlertthon W, Jo H, et al. Lung Adenocarcinoma Just Desserts: An Expanding Pie of Activating Oncogenes or a Layer Cake of Integrated Alterations. _bioRxiv_ 2025. https://doi.org/10.1101/2025.09.19.677365
