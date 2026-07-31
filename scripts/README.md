
## Dependencies

### `pPs_analysis.R`

```r
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(car)
library(vegan)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
```
Required raw data files are located in the [`../Input/pPs_Raw_data/`](../Input/pPs_Raw_data/) directory.



### `UKA_heatmap.R`

```r
library(dplyr)
library(readr)
library(pheatmap)
```
Required raw data files are located in the [`../Input/UKA_Raw_data/`](../Input/UKA_Raw_data/) directory.
### `PamDx_GTExanalysis.Rmd`
```r
library(Seurat)
library(BPCells)
library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(purrr)
library(tibble)
library(pheatmap)
library(grid)
library(data.table)
library(gridExtra)
```
Required raw data files are located in the [`../Input/GTEx_analysis/`](../Input/GTEx_analysis/) directory.

---


## Credits

The script `PamDx_GTExanalysis.Rmd` was written by [Elibrouwer](https://github.com/elibrouwer) and is included here with permission. 
