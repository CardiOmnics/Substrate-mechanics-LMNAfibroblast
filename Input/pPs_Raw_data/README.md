# Input files required for `pPs_analysis_script.R`

## Explanation
library(dplyr)    
library(pheatmap) 
library(readxl)   
library(ggplot2)  
library(readr)   

## Metadata notes
| Label | Description | Name in raw-data |
|---|---|---|
| hfLMNA<sup>WT/WT</sup> | Control human fibroblasts (wild-type *LMNA*) | WT | 
| hfLMNA<sup>R377L/WT</sup> | Cardiolaminopathic human fibroblasts (*LMNA* R377L variant) | LMNA |
| Rigid-unstrained | Conventional rigid substrate (plastic), no mechanical stimulation | Hard |
| Soft-unstrained | Soft substrate, no mechanical stimulation | Soft |
| Soft-strained | Soft substrate with cyclic stretch applied (2.5D) | Stretch |
