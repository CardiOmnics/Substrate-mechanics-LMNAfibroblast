
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
Required raw data files are located in the [`../Input/pPs_Raw_data/`](../Input/pPs_Raw_data/) directory.



### `UKA_heatmap_scripts.R`

```r
library(dplyr)    
library(pheatmap) 
library(readxl)   
library(ggplot2)  
library(readr)    
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

---


## Credits

The script `PamDx_GTExanalysis.Rmd` was written by [Elibrouwer](https://github.com/elibrouwer) and is included here with permission. 
