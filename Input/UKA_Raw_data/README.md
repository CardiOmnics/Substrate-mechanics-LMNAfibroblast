# Input required for `UKA_heatmap_scripts.R`


## Dependencies

```r
library(dplyr)    
library(pheatmap) 
library(readxl)   
library(ggplot2)  
library(readr)    
```
---

## Metadata notes

The labels are simplified/shortened in the raw-data files to improve simplicity. Please use this table as guideline through the experimental groups.

| Label | Description | Name in raw-data |
|---|---|---|
| hfLMNA<sup>WT/WT</sup> | Control human fibroblasts (wild-type *LMNA*) | WT | 
| hfLMNA<sup>R377L/WT</sup> | Cardiolaminopathic human fibroblasts (*LMNA* R377L variant) | LMNA |
| Rigid-unstrained | Conventional rigid substrate (plastic), no mechanical stimulation | Hard |
| Soft-unstrained | Soft substrate, no mechanical stimulation | Soft |
| Soft-strained | Soft substrate with cyclic stretch applied (2.5D) | Stretch |
