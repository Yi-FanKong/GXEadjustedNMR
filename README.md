# GXEadjustedNMR

`GXEadjustedNMR` is an R package for performing **G×E-adjusted nonlinear Mendelian randomization**.

The package implements a robust high-dimensional G×E-adjusted non-linear Mendelian randomization framework that first screens SNP-specific gene-environment interaction effects using `fastG×E`, then constructs a G×E-adjusted exposure and genetic instruments, and finally performs non-linear Mendelian randomization using doubly-ranked stratification.

---

## Overview

The proposed framework follows the algorithm below.

### Training set

1. For each exposure X, use `fastG×E` to screen and estimate SNP-specific G×E effects and main genetic effects.

2. Select significant G×E effects and optionally use `mmSuSiE` to prioritize causal environmental modifiers as set S.

### Analysis set

3. Construct the G×E-adjusted exposure:

$$
X^{\prime} = X - \sum_{(j,k)\in S} G_j E_k \hat{\beta}^{G \times E}_{jk}
$$

4. Construct the main-effect PRS:

$$
V^{\prime} = \sum_{j=1}^{m} G_j \hat{\beta}^{\mathrm{main}}_j
$$

5. Use X' and V' for doubly-ranked stratification.

6. Within each stratum, compute localized average causal effects, LACE:

$$
\hat{\beta}_s = \frac{\hat{\beta}^{Y}_s}{\hat{\beta}^{X}_s}
$$

7. Fit non-linear dose-response curves using fractional polynomial regression methods.

---

## Main Features

`GXEadjustedNMR` provides a set of utilities for implementing GxE-adjusted nonlinear Mendelian randomization analyses, including:

- Running `fastGxE` and `mmSuSiE` from R via `run_fastgxe()`.

> **Note:** When genotype datasets are very large, running `fastGxE` and `mmSuSiE` directly from R may fail due to excessive memory requirements. In such cases, we recommend using the original standalone implementations of `fastGxE` and `mmSuSiE`. See <https://github.com/chaoning/mmsusie> for more details.

- Selecting significant main-effect SNPs using `select_main_snps()`.
- Selecting and LD-clumping GxE SNPs using `clump_gxe_snps()`.
- Extracting genotype data using PLINK via `extract_genotypes_plink()`.
- Computing GxE adjustment scores using `compute_gxe_adjustment()`.
- Constructing GxE-adjusted exposures using `make_adjusted_exposure()`.
- Running doubly-ranked stratification using `create_nlmr_summary()`.
- Visualizing LACE estimates and heterogeneity using `plot_lace()` and `plot_heterogeneity()`.


---

## Installation

Install from GitHub:

```r
remotes::install_github("yourname/G×E-adjusted-NMR")
```

Or install from a local source directory:

```r
devtools::install("/path/to/G×EadjustedNMR")
```

Load the package:

```r
library(G×EadjustedNMR)
```

---

## External Dependencies

The full pipeline requires several external tools.

### Required command-line tools

- [`fastG×E`](https://github.com/)
- `PLINK`
- Python, if using `mmSuSiE`

### Optional Python package

- `mmsusie`

### Required R packages

The package depends on or works with:

```r
data.table
dplyr
tidyr
ggplot2
TwoSampleMR
OneSampleMR
SUMnlmr
metafor
```
