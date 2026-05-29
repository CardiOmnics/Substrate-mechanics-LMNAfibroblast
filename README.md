# Substrate mechanics modulate peptide phosphorylation and inferred kinase activity in LMNA mutant fibroblasts


## Abstract

Cells convert physical cues into biochemical signals via mechanosensitive pathways, primarily using reversible protein phosphorylation as a regulatory switch. Cardiolaminopathy caused by *LMNA* gene variants is an example of a disease with disrupted mechanosensing. Despite that, cardiolaminopathy *in vitro* modelling is still predominantly performed using rigid supraphysiological plastic and glass substrates that poorly mimic physical properties of native tissue. Therefore, we investigated how different substrates and mechanical culture conditions, here jointly referred to as substrate mechanics, regulate mechanosensitive protein phosphorylation. We profiled serine/threonine kinase (STK) activity of control and cardiolaminopathic human fibroblasts cultured across a gradient of conditions from rigid-unstrained (2D) to soft-unstrained (2D) towards an added stretch defined as soft-strained (2.5D). Principal component analysis (PCA) on the same cell lines but varying culture conditions demonstrated that substrate mechanics was the dominant factor of variability (74.4%). Following the gradient from rigid-unstrained to soft-unstrained to soft-strained conditions, we show an increase in phosphopeptide (pPs) intensity levels, improved replicate consistency, and increased inferred STK activity. We identified 20 fibroblast-expressed pPs that are substrate mechanics-responsive across the tested conditions. Furthermore, we identified unique mechanosensitive targets hidden on traditional plasticware with 55 inferred kinases in hfLMNA<sup>WT/WT</sup> and an additional 34 kinases (89 in total) in hfLMNA<sup>R377L/WT</sup> showing mechanosensitive activity when switching from a rigid to a soft substrate mechanics environment. We conclude that a biomimetic mechanical context is a key prerequisite for detecting relevant molecular drivers of cardiolaminopathy, providing an essential methodological foundation to identify pathogenic signalling pathways required to develop effective therapies.

---

## Repository Contents

| File | Description |
|---|---|
| `pPs_analysis_script.R` | Phosphopeptide (pPs) analysis: data loading and merging, heatmap, PCA, and dot plots |
| `UKA_heatmap_scripts.R` | Upstream Kinase Activity (UKA) analysis: kinase activity heatmaps for substrate and genotype comparisons |

---

## Analysis Overview

### `pPs_analysis_script.R`

This script performs the following steps:

1. **Data loading & merging** — Loads raw kinome array data from four Excel files (Runs 1–2, 3–4, and 5–6 for LMNA and WT) and merges them into a single matrix, retaining phosphopeptides present in ≥70% of samples.
2. **Heatmap** — Generates a z-score-scaled hierarchical clustering heatmap using `ComplexHeatmap`, annotated by substrate mechanics, genotype, and run ID.
3. **PCA** — Performs principal component analysis on the phosphopeptide matrix (missing values imputed with column means). Produces plots for:
   - Overall genotype separation (PC1 vs PC2)
   - Substrate mechanics separation (PC1 vs PC2)
   - Batch effects (PC2 vs PC3)
   - Variance explained bar chart (top 10 PCs)
   - Per-substrate genotype comparisons (Rigid-unstrained, Soft-unstrained, Soft-strained)
4. **Dot plots** — Visualises individual pPs expression per condition with t-tests (Bonferroni-corrected). pPs are ranked by significance; the top 25 most significant are plotted by default.

### `UKA_heatmap_scripts.R`

This script performs the following steps:

1. **Substrate mechanics comparison** — Loads UKA output CSVs comparing Hard vs Soft and Stretch vs Soft conditions, filters by contrast group, and generates a heatmap of mean kinase statistic scores (thresholded by median final score < 1.3) across substrate mechanics conditions for both genotypes.
2. **Genotype comparison** — Loads UKA output CSVs comparing hfLMNA<sup>R377L</sup> vs hfLMNA<sup>WT</sup> per substrate condition (Rigid-unstrained, Soft-unstrained, Soft-strained) and generates a corresponding kinase activity heatmap.

---

## Dependencies

### `pPs_analysis_script.R`

```r
library(readxl)         # Reading .xls input files
library(dplyr)          # Data wrangling
library(tidyr)          # pivot_longer(), unite()
library(tibble)         # column_to_rownames(), rownames_to_column()
library(purrr)          # reduce() for merging data list
library(ggplot2)        # Plotting
library(ggpubr)         # stat_compare_means(), compare_means()
library(ComplexHeatmap) # Heatmap generation
library(circlize)       # colorRamp2() for heatmap color scale
library(patchwork)      # Combining ggplot panels with | and /
```

### `UKA_heatmap_scripts.R`

```r
library(dplyr)    # Data wrangling
library(pheatmap) # Heatmap generation
library(readxl)   # Loaded but data read via readr
library(ggplot2)  # Loaded for potential downstream plotting
library(readr)    # read_csv() for UKA output files
```

---

## Input Data

Raw data files are not included in this repository. The scripts expect the following files (update paths in the scripts to match your local directory structure):

**`pPs_analysis_script.R`**
- `Raw data_Run1-2_Updated.xls`
- `Raw data_Run3-4_Updated.xls`
- `Raw data_Run5-6_LMNA_Updated.xls`
- `Raw data_Run5-6_WT_Updated.xls`

**`UKA_heatmap_scripts.R`**
- `1_UKA_Hard vs Soft Run1-2.csv`
- `2_UKA_ LMNA vs WT run1-2.csv`
- `3_UKA_Stretch vs Soft Run3-4.csv`
- `4_UKA_LMNA vs WT run3-4.csv`

---

## Experimental Groups

| Label | Description |
|---|---|
| `hfLMNA_WT` | Control human fibroblasts (wild-type *LMNA*) |
| `hfLMNA_R377L` | Cardiolaminopathic human fibroblasts (*LMNA* R377L variant) |
| Rigid-unstrained | Conventional rigid substrate (plastic/glass), no mechanical stimulation |
| Soft-unstrained | Soft hydrogel substrate, no mechanical stimulation |
| Soft-strained | Soft hydrogel substrate with cyclic stretch applied (2.5D) |
