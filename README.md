# Substrate mechanics modulate peptide phosphorylation and inferred kinase activity in LMNA mutant fibroblasts


## Abstract

Cells convert physical cues into biochemical signals via mechanosensitive pathways, primarily using reversible protein phosphorylation as a regulatory switch. Cardiolaminopathy caused by *LMNA* gene variants is an example of a disease with disrupted mechanosensing. Despite that, cardiolaminopathy *in vitro* modelling is still predominantly performed using rigid supraphysiological plastic and glass substrates that poorly mimic physical properties of native tissue. Therefore, we investigated how different substrates and mechanical culture conditions, here jointly referred to as substrate mechanics, regulate mechanosensitive protein phosphorylation. We profiled serine/threonine kinase (STK) activity of control and cardiolaminopathic human fibroblasts cultured across a gradient of conditions from rigid-unstrained (2D) to soft-unstrained (2D) towards an added stretch defined as soft-strained (2.5D). Principal component analysis (PCA) on the same cell lines but varying culture conditions demonstrated that substrate mechanics was the dominant factor of variability (74.4%). Following the gradient from rigid-unstrained to soft-unstrained to soft-strained conditions, we show an increase in phosphopeptide (pPs) intensity levels, improved replicate consistency, and increased inferred STK activity. We identified 20 fibroblast-expressed pPs that are substrate mechanics-responsive across the tested conditions. Furthermore, we identified unique mechanosensitive targets hidden on traditional plasticware with 55 inferred kinases in hfLMNA<sup>WT/WT</sup> and an additional 34 kinases (89 in total) in hfLMNA<sup>R377L/WT</sup> showing mechanosensitive activity when switching from a rigid to a soft substrate mechanics environment. We conclude that a biomimetic mechanical context is a key prerequisite for detecting relevant molecular drivers of cardiolaminopathy, providing an essential methodological foundation to identify pathogenic signalling pathways required to develop effective therapies.

---

## Repository Contents

| File | Description |
|---|---|
| `pPs_analysis_script.R` | Phosphopeptide (pPs) analysis: data loading and merging, heatmap, PCA, and dot plots |
| `UKA_heatmap_scripts.R` | Upstream Kinase Activity (UKA) analysis: kinase activity heatmaps for substrate and genotype comparisons |
| `PamDx_GTExanalysis.Rmd` | GTEx analysis: Validation of pPs expression in human fibroblast (skin and cultured) |

---

## Analysis Overview

### `pPs_analysis_script.R`

This script performs the following steps:

1. **Data loading & merging** — Loads raw kinome array data from four Excel files (Runs 1–2, 3–4, and 5–6 for LMNA and WT) and merges them into a single matrix, retaining phosphopeptides present in ≥90% of samples.
2. **Heatmap** — Generates a z-score-scaled hierarchical clustering heatmap using `ComplexHeatmap`, annotated by substrate mechanics, genotype, and run ID.
3. **PCA** — Performs principal component analysis on the phosphopeptide matrix (missing values imputed with column means). Produces plots for:
   - Overall genotype separation (_PC1 vs PC2_)
   - Substrate mechanics separation (_PC1 vs PC2_)
   - Batch effects (_PC2 vs PC3_)
   - Variance explained bar chart (top 10 PCs)
   - Per-substrate genotype comparisons (Rigid-unstrained, Soft-unstrained, Soft-strained)
4. **Dot plots** — Visualises individual pPs expression per condition with t-tests (Bonferroni-corrected). pPs are ranked by significance; the top 25 most significant are plotted by default.

### `UKA_heatmap_scripts.R`

This script performs the following steps:

1. **Substrate mechanics comparison** — Loads UKA output CSVs comparing Hard vs Soft and Stretch vs Soft conditions, filters by contrast group, and generates a heatmap of mean kinase statistic scores (thresholded by median final score < 1.3) across substrate mechanics conditions for both genotypes.
2. **Genotype comparison** — Loads UKA output CSVs comparing hfLMNA<sup>R377L/WT</sup> vs hfLMNA<sup>WT/WT</sup> per substrate condition (Rigid-unstrained, Soft-unstrained, Soft-strained) and generates a corresponding kinase activity heatmap.

---

## Dependencies

### `pPs_analysis_script.R`

```r
library(readxl)        
library(dplyr)          
library(tidyr)         
library(tibble)         
library(purrr)          
library(ggplot2)        
library(ggpubr)         
library(ComplexHeatmap) 
library(circlize)       
library(patchwork)      
```

### `UKA_heatmap_scripts.R`

```r
library(dplyr)    
library(pheatmap) 
library(readxl)   
library(ggplot2)  
library(readr)    
```

---

## Input Data

The scripts expect the following files (update paths in the scripts to match your local directory):

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

The labels are simplified/shortened in the raw-data files to improve simplicity. Please use this table as guideline through the experimental groups.

| Label | Description | Name in raw-data |
|---|---|---|
| hfLMNA<sup>WT/WT</sup> | Control human fibroblasts (wild-type *LMNA*) | WT | 
| hfLMNA<sup>R377L/WT</sup> | Cardiolaminopathic human fibroblasts (*LMNA* R377L variant) | LMNA |
| Rigid-unstrained | Conventional rigid substrate (plastic), no mechanical stimulation | Hard |
| Soft-unstrained | Soft substrate, no mechanical stimulation | Soft |
| Soft-strained | Soft substrate with cyclic stretch applied (2.5D) | Stretch |


## Credits

The script PamDx_GTExanalysis.Rmd was written by @ and is included here with permission. 

