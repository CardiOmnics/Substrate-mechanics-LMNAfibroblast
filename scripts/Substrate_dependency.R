# ==============================================================================
# hfLMNA (R377L) kinome / phosphopeptide (PamStation) analysis
# ==============================================================================

## ---- 0. Setup ---------------------------------------------------------------
## IF NOT INSTALLED:

# install.packages(c("readxl","dplyr","tidyr","tibble","purrr","ggplot2",
#                     "ggpubr","ggsignif","ComplexHeatmap","circlize","patchwork",
#                     "car","vegan","lme4","lmerTest","emmeans","performance"))


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

surface_order <- c("Rigid-unstrained", "Soft-unstrained", "Soft-strained")

genotype_colors <- c("hfLMNA_WT" = "chartreuse4", "hfLMNA_R377L" = "orangered2")
surface_shapes  <- c(`Rigid-unstrained` = 16, `Soft-unstrained` = 17, `Soft-strained` = 15)
surface_colors  <- c("Rigid-unstrained" = "#ee7faaff",
                     "Soft-unstrained" = "#b4a7d6ff",
                     "Soft-strained"   = "#6d9eebff")
condition_surface_colors <- c(
  "hfLMNA_WT Rigid-unstrained"     = "chartreuse4",
  "hfLMNA_WT Soft-unstrained"      = "chartreuse3",
  "hfLMNA_WT Soft-strained"        = "chartreuse2",
  "hfLMNA_R377L Rigid-unstrained"  = "orangered3",
  "hfLMNA_R377L Soft-unstrained"   = "orangered2",
  "hfLMNA_R377L Soft-strained"     = "orangered1"
)

common_theme <- theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = -10),
        legend.spacing.y = unit(0.5, "lines"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))

# ---- File paths (edit to match your local data layout) --------------
data_dir <- "YOUR/DATA/REPOSITORY/"
# ----
f_run1_2      <- file.path(data_dir, "Raw data_Run1-2_Updated.xls")
f_run3_4      <- file.path(data_dir, "Raw data_Run3-4_Updated.xls")
f_run5_6_lmna <- file.path(data_dir, "Raw data_Run5-6_LMNA_Updated.xls")
f_run5_6_wt   <- file.path(data_dir, "Raw data_Run5-6_WT_Updated.xls")


## ---- 1. Import & merge raw kinome data ---------------------------------------

Raw_1_2      <- read_excel(f_run1_2, skip = 4)
Raw_3_4      <- read_excel(f_run3_4, skip = 4)
Raw_5_6_LMNA <- read_excel(f_run5_6_lmna, skip = 4)
Raw_5_6_WT   <- read_excel(f_run5_6_wt, skip = 4)

Raw_1_2 <- Raw_1_2 %>%
  rename_with(~ paste0(., "_RUN1_2"), -c(Description, ID, Sequence, `...4`))
Raw_3_4 <- Raw_3_4 %>%
  rename_with(~ paste0(., "_RUN3_4"), -c(Description, ID, Sequence, `...4`))
Raw_5_6_LMNA <- Raw_5_6_LMNA %>%
  rename_with(~ paste0(., "_RUN5_6"), -c(Description, ID, Sequence, `...4`))
Raw_5_6_WT <- Raw_5_6_WT %>%
  rename_with(~ paste0(., "_RUN5_6"), -c(Description, ID, Sequence, `...4`))
Raw_5_6_LMNA_HardCtrl <- Raw_5_6_LMNA %>% dplyr::select(1:3, contains("Hard Ctrl"))
Raw_5_6_WT_HardCtrl   <- Raw_5_6_WT   %>% dplyr::select(1:3, contains("Hard Ctrl"))

raw_combined_data <- purrr::reduce(
  list(Raw_1_2, Raw_3_4, Raw_5_6_LMNA_HardCtrl, Raw_5_6_WT_HardCtrl),
  full_join,
  by = colnames(Raw_1_2)[1:3]
)

heatmap_df <- raw_combined_data %>%
  unite("ID", 2, remove = TRUE) %>%
  dplyr::select(-contains("...4"), -contains("4.y"))

heatmap_matrix <- heatmap_df %>%
  column_to_rownames("ID") %>%
  dplyr::select(-1, -2) %>%
  as.matrix()

keep_rows <- rowSums(!is.na(heatmap_matrix)) / ncol(heatmap_matrix) >= 0.9
heatmap_matrix <- heatmap_matrix[keep_rows, ]


## ---- 2. Sample annotation -----------------------------------------------------

Sub_annotation_df <- data.frame(
  Substrate_mechanics = factor(c(
    "Soft-unstrained", "Rigid-unstrained", "Soft-unstrained", "Rigid-unstrained", "Soft-unstrained", "Rigid-unstrained", "Soft-unstrained", "Rigid-unstrained",
    "Soft-unstrained", "Soft-unstrained", "Rigid-unstrained", "Soft-unstrained", "Soft-strained", "Soft-strained", "Soft-unstrained", "Soft-strained",
    "Soft-unstrained", "Soft-strained", "Soft-unstrained", "Soft-strained", "Soft-strained", "Soft-unstrained", "Soft-strained", "Soft-strained",
    "Soft-unstrained", "Rigid-unstrained", "Rigid-unstrained", "Rigid-unstrained", "Rigid-unstrained", "Rigid-unstrained", "Rigid-unstrained", "Rigid-unstrained", "Soft-unstrained",
    "Soft-unstrained", "Soft-unstrained", "Soft-unstrained", "Rigid-unstrained", "Rigid-unstrained", "Rigid-unstrained",
    "Rigid-unstrained"
  ), levels = surface_order),
  Genotype = factor(c(
    "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_R377L",
    "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L",
    "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT",
    "hfLMNA_WT", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L",
    "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_WT",
    "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_R377L", "hfLMNA_WT", "hfLMNA_R377L",
    "hfLMNA_WT", "hfLMNA_WT", "hfLMNA_R377L", "hfLMNA_R377L", "hfLMNA_WT"
  )),
  RunID = factor(c(
    "Run1_2", "Run1_2", "Run1_2",
    "Run1_2", "Run1_2", "Run1_2",
    "Run1_2", "Run1_2", "Run1_2",
    "Run1_2", "Run1_2", "Run3_4",
    "Run3_4", "Run3_4", "Run3_4",
    "Run3_4", "Run3_4", "Run3_4",
    "Run3_4", "Run3_4", "Run3_4",
    "Run3_4", "Run3_4", "Run3_4",
    "Run3_4", "Run5_6", "Run5_6",
    "Run5_6", "Run5_6", "Run5_6",
    "Run5_6", "Run1_2", "Run3_4",
    "Run1_2", "Run1_2", "Run3_4",
    "Run5_6", "Run5_6", "Run1_2",
    "Run1_2"
  )),
  row.names = c(
    "27. WT Soft Ctrl 3_RUN1_2", "1. WT Hard Ctrl 1_RUN1_2", "25. WT Soft Ctrl 1_RUN1_2",
    "2. WT Hard Ctrl 2_RUN1_2", "26. WT Soft Ctrl 2_RUN1_2", "15. LMNA Hard Ctrl 3_RUN1_2",
    "39. LMNA Soft Ctrl 3_RUN1_2", "13. LMNA Hard Ctrl 1_RUN1_2", "37. LMNA Soft Ctrl 1_RUN1_2",
    "38. LMNA Soft Ctrl 2_RUN1_2", "14. LMNA Hard Ctrl 2_RUN1_2", "27. WT Soft Ctrl 3_RUN3_4",
    "51. WT Stretch Ctrl 3_RUN3_4", "52. WT Stretch Ctrl 4_RUN3_4", "25. WT Soft Ctrl 1_RUN3_4",
    "49. WT Stretch Ctrl 1_RUN3_4", "26. WT Soft Ctrl 2_RUN3_4", "50. WT Stretch Ctrl 2_RUN3_4",
    "39. LMNA Soft Ctrl 3_RUN3_4", "63. LMNA Stretch Ctrl 3_RUN3_4", "64. LMNA Stretch Ctrl 4_RUN3_4",
    "37. LMNA Soft Ctrl 1_RUN3_4", "61. LMNA Stretch Ctrl 1_RUN3_4", "61. LMNA Stretch Ctrl 2_RUN3_4",
    "38. LMNA Soft Ctrl 1_RUN3_4", "13. LMNA Hard Ctrl 1_RUN5_6", "14. LMNA Hard Ctrl 2_RUN5_6",
    "15. LMNA Hard Ctrl 3_RUN5_6", "1. WT Hard Ctrl 1_RUN5_6", "2. WT Hard Ctrl 2_RUN5_6",
    "3. WT Hard Ctrl 3_RUN5_6", "4. WT Hard Ctrl 4_RUN1_2", "40. LMNA Soft Ctrl 4_RUN3_4",
    "28. WT Soft Ctrl 4_RUN1_2", "40. LMNA Soft Ctrl 4_RUN1_2", "28. WT Soft Ctrl 4_RUN3_4",
    "4. WT Hard Ctrl 4_RUN5_6", "16. LMNA Hard Ctrl 4_RUN5_6", "16. LMNA Hard Ctrl 4_RUN1_2",
    "3. WT Hard Ctrl 3_RUN1_2"
  )
)

Sub_ann_colors <- list(
  Substrate_mechanics = surface_colors,
  Genotype = genotype_colors,
  RunID = c("Run1_2" = "grey85", "Run3_4" = "grey65", "Run5_6" = "grey40")
)

new_order <- c(
  "1. WT Hard Ctrl 1_RUN1_2", "2. WT Hard Ctrl 2_RUN1_2", "3. WT Hard Ctrl 3_RUN1_2", "4. WT Hard Ctrl 4_RUN1_2",
  "1. WT Hard Ctrl 1_RUN5_6", "2. WT Hard Ctrl 2_RUN5_6", "3. WT Hard Ctrl 3_RUN5_6", "4. WT Hard Ctrl 4_RUN5_6",
  
  "25. WT Soft Ctrl 1_RUN1_2", "26. WT Soft Ctrl 2_RUN1_2", "27. WT Soft Ctrl 3_RUN1_2", "28. WT Soft Ctrl 4_RUN1_2",
  "25. WT Soft Ctrl 1_RUN3_4", "26. WT Soft Ctrl 2_RUN3_4", "27. WT Soft Ctrl 3_RUN3_4", "28. WT Soft Ctrl 4_RUN3_4",
  
  "49. WT Stretch Ctrl 1_RUN3_4", "50. WT Stretch Ctrl 2_RUN3_4", "51. WT Stretch Ctrl 3_RUN3_4", "52. WT Stretch Ctrl 4_RUN3_4",
  
  "13. LMNA Hard Ctrl 1_RUN1_2", "14. LMNA Hard Ctrl 2_RUN1_2", "15. LMNA Hard Ctrl 3_RUN1_2", "16. LMNA Hard Ctrl 4_RUN1_2",
  "13. LMNA Hard Ctrl 1_RUN5_6", "14. LMNA Hard Ctrl 2_RUN5_6", "15. LMNA Hard Ctrl 3_RUN5_6", "16. LMNA Hard Ctrl 4_RUN5_6",
  
  "37. LMNA Soft Ctrl 1_RUN1_2", "38. LMNA Soft Ctrl 2_RUN1_2", "39. LMNA Soft Ctrl 3_RUN1_2", "40. LMNA Soft Ctrl 4_RUN1_2",
  "37. LMNA Soft Ctrl 1_RUN3_4", "38. LMNA Soft Ctrl 1_RUN3_4", "39. LMNA Soft Ctrl 3_RUN3_4", "40. LMNA Soft Ctrl 4_RUN3_4",
  
  "61. LMNA Stretch Ctrl 1_RUN3_4", "61. LMNA Stretch Ctrl 2_RUN3_4", "63. LMNA Stretch Ctrl 3_RUN3_4", "64. LMNA Stretch Ctrl 4_RUN3_4"
)
heatmap_matrix <- heatmap_matrix[, new_order]
design_df <- Sub_annotation_df[colnames(heatmap_matrix), ]


## ---- 3. Heatmap ----------------------------------------------------------------

col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

ComplexHeatmap::pheatmap(
  t(heatmap_matrix),
  scale = "column",
  col = col_fun,
  annotation_row = Sub_annotation_df,
  annotation_colors = Sub_ann_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  column_split = 7,
  na_col = "grey80",
  gaps_row = c(8, 16, 20, 28, 36),
  show_rownames = FALSE,
  heatmap_legend_param = list(title = "Z-score",
                              title_position = "topcenter",
                              legend_direction = "horizontal")
) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom")


## ---- 4. PCA -----------------------------------------------------------------
pca_matrix <- apply(
  heatmap_matrix, 2,
  function(x) {
    x[!is.finite(x)] <- NA
    replace(x, is.na(x), mean(x, na.rm = TRUE))
  }
)

pca <- prcomp(t(pca_matrix), scale. = TRUE, center = TRUE)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3]
)
pca_df$Condition <- ifelse(grepl("LMNA", pca_df$Sample), "hfLMNA_R377L", "hfLMNA_WT")
pca_df$Surface <- case_when(
  grepl("Hard", pca_df$Sample, ignore.case = TRUE)    ~ "Rigid-unstrained",
  grepl("Soft", pca_df$Sample, ignore.case = TRUE)    ~ "Soft-unstrained",
  grepl("Stretch", pca_df$Sample, ignore.case = TRUE) ~ "Soft-strained"
)
pca_df$Batch <- case_when(
  grepl("RUN1_2", pca_df$Sample) ~ "R1-2",
  grepl("RUN3_4", pca_df$Sample) ~ "R3-4",
  grepl("RUN5_6", pca_df$Sample) ~ "R5-6",
  TRUE ~ "Other"
)
pca_df$Condition <- factor(pca_df$Condition, levels = c("hfLMNA_WT", "hfLMNA_R377L"))
pca_df$Surface   <- factor(pca_df$Surface, levels = surface_order)

pc1_lab <- paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")
pc3_lab <- paste0("PC3 (", round(summary(pca)$importance[2, 3] * 100, 1), "%)")

### 4.1 Genotype (LMNA vs WT), all substrates combined ------------------------

centroids_p1 <- pca_df %>%
  group_by(Condition) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
d_p1 <- as.numeric(dist(centroids_p1[, c("PC1", "PC2")]))

p1 <- ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Surface)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Condition), linewidth = 0.9, type = "t") +
  geom_segment(data = centroids_p1,
               aes(x = PC1[1], y = PC2[1], xend = PC1[2], yend = PC2[2]),
               inherit.aes = FALSE, color = "black") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3,
           label = paste0("italic(d) ==", round(d_p1, 2)), parse = TRUE) +
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
  scale_shape_manual(values = surface_shapes, breaks = surface_order) +
  guides(shape = "none") +
  labs(title = "", color = "", x = pc1_lab, y = pc2_lab)

### 4.2 All substrate mechanics ------------------------------------

centroids_p2 <- pca_df %>%
  group_by(Surface) %>%
  summarise(PC1 = mean(PC1, na.rm = TRUE), PC2 = mean(PC2, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()
rownames(centroids_p2) <- centroids_p2$Surface
dist_matrix <- as.matrix(dist(centroids_p2[, c("PC1", "PC2")]))
d_RS <- dist_matrix["Rigid-unstrained", "Soft-unstrained"]
d_RT <- dist_matrix["Rigid-unstrained", "Soft-strained"]
d_ST <- dist_matrix["Soft-unstrained", "Soft-strained"]

p2 <- ggplot(pca_df, aes(PC1, PC2, color = Surface, shape = Surface)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Surface), linewidth = 0.9, type = "t", level = 0.95) +
  geom_segment(aes(x = centroids_p2$PC1[1], y = centroids_p2$PC2[1],
                   xend = centroids_p2$PC1[2], yend = centroids_p2$PC2[2]),
               inherit.aes = FALSE, color = "black", linewidth = 0.7) +
  geom_segment(aes(x = centroids_p2$PC1[1], y = centroids_p2$PC2[1],
                   xend = centroids_p2$PC1[3], yend = centroids_p2$PC2[3]),
               inherit.aes = FALSE, color = "black", linewidth = 0.7) +
  geom_segment(aes(x = centroids_p2$PC1[2], y = centroids_p2$PC2[2],
                   xend = centroids_p2$PC1[3], yend = centroids_p2$PC2[3]),
               inherit.aes = FALSE, color = "black", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.1, size = 3,
           label = paste0(
             "italic(d)(RU-SS)==", round(d_RT, 2), "*','~",
             "italic(d)(RU-SU)==", round(d_RS, 2), "*','~",
             "italic(d)(SU-SS)==", round(d_ST, 2)
           ), parse = TRUE) +
  theme_classic() +
  scale_shape_manual(values = surface_shapes, breaks = surface_order) +
  scale_color_manual(values = surface_colors, breaks = surface_order) +
  guides(shape = "none") +
  labs(title = "", color = "", x = pc1_lab, y = pc2_lab)

### 4.3 Batch effect (PC2 vs PC3) ----------------------------------------

centroids_p3 <- pca_df %>%
  group_by(Batch) %>%
  summarise(PC2 = mean(PC2, na.rm = TRUE), PC3 = mean(PC3, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()
rownames(centroids_p3) <- centroids_p3$Batch
dist_matrix2 <- as.matrix(dist(centroids_p3[, c("PC2", "PC3")]))
d_14 <- dist_matrix2["R1-2", "R3-4"]
d_16 <- dist_matrix2["R1-2", "R5-6"]
d_36 <- dist_matrix2["R3-4", "R5-6"]

p3 <- ggplot(pca_df, aes(PC2, PC3, color = Batch, shape = Surface)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Batch), linewidth = 0.9, type = "t", level = 0.95) +
  geom_segment(aes(x = centroids_p3$PC2[1], y = centroids_p3$PC3[1],
                   xend = centroids_p3$PC2[2], yend = centroids_p3$PC3[2]),
               inherit.aes = FALSE, color = "black", linewidth = 0.7) +
  geom_segment(aes(x = centroids_p3$PC2[1], y = centroids_p3$PC3[1],
                   xend = centroids_p3$PC2[3], yend = centroids_p3$PC3[3]),
               inherit.aes = FALSE, color = "black", linewidth = 0.7, linetype = "dashed") +
  geom_segment(aes(x = centroids_p3$PC2[2], y = centroids_p3$PC3[2],
                   xend = centroids_p3$PC2[3], yend = centroids_p3$PC3[3]),
               inherit.aes = FALSE, color = "black", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.1, size = 3,
           label = paste0(
             "italic(d)(R12~R56)==", round(d_16, 2), "*','~",
             "italic(d)(R12~R34)==", round(d_14, 2), "*','~",
             "italic(d)(R34~R56)==", round(d_36, 2)
           ), parse = TRUE) +
  theme_classic() +
  scale_color_manual(values = c(`R1-2` = "grey85", `R3-4` = "grey65", `R5-6` = "grey40", `R7-8` = "grey15")) +
  scale_shape_manual(values = surface_shapes, breaks = surface_order) +
  guides(shape = "none") +
  labs(title = "", color = "", x = pc2_lab, y = pc3_lab) +
  theme(plot.margin = margin(t = 1, r = 2, b = 2, l = 6))

### 4.4 Variance explained ---------------------------------------------

pca_var <- summary(pca)$importance[2, ] * 100
pca_var_df <- data.frame(PC = paste0("PC", seq_along(pca_var)), Variance = as.numeric(pca_var))

p4 <- ggplot(pca_var_df[1:10, ], aes(x = factor(PC, levels = PC[1:10]), y = Variance)) +
  geom_col(fill = "black", alpha = 0.8, width = 0.9) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")),
            vjust = -0.3, hjust = 0.5, size = 2.5, fontface = "bold", check_overlap = TRUE) +
  labs(title = "", x = "", y = "% Variance explained") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 1, r = 2, b = 2, l = 6)
  )

p1_th <- p1 + common_theme
p2_th <- p2 + common_theme
p3_th <- p3 + common_theme
p4_th <- p4 + common_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

### 4.5 Genotype separation in each substrate mechanic -------------------

plot_genotype_within_surface <- function(surface_level) {
  pdat <- pca_df[pca_df$Surface == surface_level, ]
  pdat$Condition <- factor(pdat$Condition, levels = c("hfLMNA_WT", "hfLMNA_R377L"))
  
  centroids <- pdat %>%
    group_by(Condition) %>%
    summarise(PC1 = mean(PC1, na.rm = TRUE), PC2 = mean(PC2, na.rm = TRUE), .groups = "drop")
  d_euclid <- as.numeric(dist(centroids[, c("PC1", "PC2")], method = "euclidean"))
  
  ggplot(pdat, aes(PC1, PC2, color = Condition, shape = Surface)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Condition), linewidth = 0.9, type = "t") +
    geom_segment(data = centroids,
                 aes(x = PC1[1], y = PC2[1], xend = PC1[2], yend = PC2[2]),
                 inherit.aes = FALSE, color = "black", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3, color = "black",
             label = paste0("italic(d) ==", round(d_euclid, 2)), parse = TRUE) +
    theme_classic() +
    scale_color_manual(values = genotype_colors) +
    scale_shape_manual(values = surface_shapes, breaks = surface_order) +
    guides(shape = "none") +
    labs(title = "", color = "", x = pc1_lab, y = pc2_lab) +
    common_theme
}

p1_RIGID   <- plot_genotype_within_surface("Rigid-unstrained")
p1_SOFT    <- plot_genotype_within_surface("Soft-unstrained")
p1_STRETCH <- plot_genotype_within_surface("Soft-strained")

((p1_th | p2_th | p3_th / p4_th) / (p1_RIGID | p1_SOFT | p1_STRETCH))


## ---- 5. Variance partitioning & PERMANOVA -------------------------------------

pc_df <- cbind(
  as.data.frame(pca$x[, 1:5]),
  design_df[rownames(pca$x), c("Substrate_mechanics", "Genotype", "RunID")]
)

pc_variance_partition <- purrr::map_dfr(paste0("PC", 1:5), function(pc) {
  fit <- lm(reformulate(c("Substrate_mechanics", "Genotype", "RunID"), response = pc), data = pc_df)
  tab <- car::Anova(fit, type = "II")
  data.frame(
    PC = pc, Term = rownames(tab)[1:3],
    EtaSq = tab$`Sum Sq`[1:3] / sum(tab$`Sum Sq`),
    p = tab$`Pr(>F)`[1:3]
  )
})
print(pc_variance_partition)

dist_mat <- dist(t(pca_matrix), method = "euclidean")

set.seed(1)
permanova_full <- adonis2(
  dist_mat ~ Substrate_mechanics + Genotype + RunID,
  data = design_df, permutations = 999, by = "margin"
)
print(permanova_full)
anova(betadisper(dist_mat, design_df$Substrate_mechanics))
anova(betadisper(dist_mat, design_df$Genotype))


## ---- 6. Peptide-level differential phosphorylation (dot plots) ---------------

### 6.1 Long-format data + condition/surface labels ----------------------------

dot_data <- heatmap_matrix %>%
  as.data.frame() %>%
  rownames_to_column("pPs") %>%
  pivot_longer(-pPs, names_to = "Sample", values_to = "Expression") %>%
  mutate(
    Condition = case_when(
      grepl("LMNA *Hard", Sample, ignore.case = TRUE)    ~ "hfLMNA_R377L Rigid-unstrained",
      grepl("LMNA_Hard", Sample, ignore.case = TRUE)      ~ "hfLMNA_R377L Rigid-unstrained",
      grepl("WT *Hard", Sample, ignore.case = TRUE)       ~ "hfLMNA_WT Rigid-unstrained",
      grepl("WT_Hard", Sample, ignore.case = TRUE)        ~ "hfLMNA_WT Rigid-unstrained",
      grepl("LMNA *Soft", Sample, ignore.case = TRUE)     ~ "hfLMNA_R377L Soft-unstrained",
      grepl("WT *Soft", Sample, ignore.case = TRUE)       ~ "hfLMNA_WT Soft-unstrained",
      grepl("LMNA *Stretch", Sample, ignore.case = TRUE)  ~ "hfLMNA_R377L Soft-strained",
      grepl("WT *Stretch", Sample, ignore.case = TRUE)    ~ "hfLMNA_WT Soft-strained",
      TRUE ~ "Other"
    ),
    CellGroup = case_when(
      grepl("LMNA", Sample, ignore.case = TRUE) ~ "hfLMNA_R377L",
      grepl("WT", Sample, ignore.case = TRUE)   ~ "hfLMNA_WT",
      TRUE ~ "Other"
    ),
    Surface = sub("^.*?\\s", "", Condition)
  ) %>%
  filter(Condition != "Other")

dot_data$Condition <- factor(dot_data$Condition, levels = names(condition_surface_colors))
dot_data$Surface   <- factor(dot_data$Surface, levels = surface_order)

# Peptides excluded from analysis (e.g. flagged as unreliable/QC failures)
exclude_pPs <- unique(c(
  "EPB42_241_253", "ACM5_494_506", "KPCB_19_31_A25S", "SCN7A_898_910",
  "CFTR_730_742", "CFTR_761_773", "ADRB2_338_350", "RYR1_4317_4329",
  "PLM_76_88", "MYPC3_268_280", "TOP2A_1463_1475", "ADDB_696_708",
  "ADDB_706_718", "NCF1_296_308", "NCF1_321_333", "GBRB2_427_439",
  "NOS3_1171_1183", "PLEK_106_118", "GYS2_1_13", "PTK6_436_448",
  "NMDZ1_890_902", "KCNA2_442_454", "MPIP1_172_184", "MPIP3_208_220",
  "PPR1A_28_40", "GSUB_61_73", "KCNA3_461_473", "RS6_228_240",
  "NEK2_172_184", "KCNA1_438_450", "GPR6_349_361", "CGHB_109_121",
  "KCNA6_504_516", "pTY3H_64_78", "TY3H_65_77", "ACM4_456_468",
  "ACM5_498_510"
))
dot_data <- dot_data %>% filter(!pPs %in% exclude_pPs)

comparisons <- list(
  c("hfLMNA_WT Rigid-unstrained", "hfLMNA_WT Soft-unstrained"),
  c("hfLMNA_WT Soft-unstrained", "hfLMNA_WT Soft-strained"),
  c("hfLMNA_R377L Rigid-unstrained", "hfLMNA_R377L Soft-unstrained"),
  c("hfLMNA_R377L Soft-unstrained", "hfLMNA_R377L Soft-strained")
)

### 6.2 Multiplicity-corrected peptide statistics -------------------------------
raw_pvals <- dot_data %>%
  group_by(pPs) %>%
  group_modify(~ compare_means(
    Expression ~ Condition, data = .x, method = "t.test",
    comparisons = comparisons, p.adjust.method = "none"
  )) %>%
  ungroup()

pval_stats_df <- raw_pvals %>%
  group_by(group1, group2) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(p.signif = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

pval_stats <- pval_stats_df %>%
  group_by(pPs) %>%
  summarise(min_pval = min(p.adj, na.rm = TRUE),
            n_sig = sum(p.adj < 0.05, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_pval, desc(n_sig)) %>%
  pull(pPs)

top25_pPs <- pval_stats[1:25]
rest_pPs  <- pval_stats[!(pval_stats %in% top25_pPs)]

### 6.3 Reusable plotting helpers ------------------------------------------------

build_signif_annotations <- function(pPs_subset) {
  pep_range <- dot_data %>%
    filter(pPs %in% pPs_subset) %>%
    group_by(pPs) %>%
    summarise(y_max = max(Expression, na.rm = TRUE),
              y_range = diff(range(Expression, na.rm = TRUE)), .groups = "drop")
  
  pval_stats_df %>%
    filter(pPs %in% pPs_subset) %>%
    mutate(pPs = factor(pPs, levels = pPs_subset)) %>%
    left_join(pep_range, by = "pPs") %>%
    group_by(pPs) %>%
    mutate(
      xmin = match(group1, levels(dot_data$Condition)),
      xmax = match(group2, levels(dot_data$Condition)),
      y_position = y_max + 0.08 * y_range + (row_number() - 1) * 0.14 * y_range
    ) %>%
    ungroup()
}

plot_pPs_panel <- function(pPs_subset, ncol = 5) {
  plot_data <- dot_data %>%
    filter(pPs %in% pPs_subset) %>%
    mutate(pPs = factor(pPs, levels = pPs_subset))
  annot_df <- build_signif_annotations(pPs_subset)
  
  ggplot(plot_data, aes(x = Condition, y = Expression, color = Condition, shape = Surface)) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 2) +
    geom_smooth(aes(group = CellGroup), method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    stat_compare_means(
      comparisons = comparisons,
      method = "t.test",
      p.adjust.method = "BH",
      label = "p.signif",
      step.increase = 0.05,
      size = 2.5
    )  +
    facet_wrap(~ pPs, ncol = ncol, scales = "free") +
    scale_color_manual(values = condition_surface_colors, guide = guide_legend(nrow = 3, byrow = FALSE)) +
    scale_shape_manual(values = surface_shapes, breaks = surface_order) +
    scale_y_continuous(breaks = function(x) unique(round(pretty(x, n = 4))),
                       expand = expansion(mult = c(0.05, 0.25))) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", axis.text.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black")) +
    guides(shape = "none") +
    labs(y = "QCLog_CmbCo intensity", x = "")
}

### 6.4 Panels: top 25 most significant peptides, and the remainder ------------

p_top25 <- plot_pPs_panel(top25_pPs, ncol = 5)
p_top25

p_rest <- plot_pPs_panel(rest_pPs, ncol = 5)
p_rest


## ---- 7. Linear mixed model: genotype x mechanical condition interaction ----

lmm_data <- dot_data %>%
  mutate(
    Genotype = factor(CellGroup, levels = c("hfLMNA_WT", "hfLMNA_R377L")),
    Substrate_mechanics = factor(Surface, levels = surface_order),  # unordered - no gradient assumption
    pPs = factor(pPs),
    Sample = factor(Sample)
  )

lmm_fit <- lmer(
  Expression ~ Substrate_mechanics * Genotype + (1 | pPs) + (1 | Sample),
  data = lmm_data
)
summary(lmm_fit)
anova(lmm_fit)                
performance::r2(lmm_fit)

emm <- emmeans(lmm_fit, ~ Substrate_mechanics | Genotype)
planned_contrasts <- contrast(emm, method = list(
  "Rigid vs Soft"   = c(-1, 1, 0),
  "Soft vs Stretch" = c(0, -1, 1)
), adjust = "BH")
summary(planned_contrasts)

### 7.1 Model-estimated means plot ----------------------------------------------

emm_df <- as.data.frame(emm) %>%
  mutate(Condition = factor(paste(Genotype, Substrate_mechanics), levels = levels(dot_data$Condition)))

contrast_df <- as.data.frame(planned_contrasts) %>%
  mutate(
    group_from = ifelse(contrast == "Rigid vs Soft", "Rigid-unstrained", "Soft-unstrained"),
    group_to   = ifelse(contrast == "Rigid vs Soft", "Soft-unstrained", "Soft-strained"),
    xmin = match(paste(Genotype, group_from), levels(dot_data$Condition)),
    xmax = match(paste(Genotype, group_to), levels(dot_data$Condition)),
    p.signif = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                         p.value < 0.05 ~ "*", TRUE ~ "ns")
  ) %>%
  arrange(Genotype) %>%
  mutate(y_position = max(dot_data$Expression, na.rm = TRUE) +
           0.04 * diff(range(dot_data$Expression, na.rm = TRUE)) +
           (row_number() - 0.5) * 0.05 * diff(range(dot_data$Expression, na.rm = TRUE)))

ggplot(dot_data, aes(x = Condition, y = Expression)) +
  geom_jitter(aes(color = Condition, shape = Surface), width = 0.2, alpha = 0.8, size = 2) +
  geom_line(data = emm_df, aes(x = Condition, y = emmean, group = Genotype),
            linewidth = 0.7, color = "black") +
  geom_errorbar(data = emm_df, aes(x = Condition, y = emmean, ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15, linewidth = 0.6, color = "black") +
  geom_point(data = emm_df, aes(x = Condition, y = emmean),
             size = 3, color = "black", fill = "black", shape = 21, stroke = 1) +
  geom_signif(data = contrast_df,
              aes(xmin = xmin, xmax = xmax, annotations = p.signif, y_position = y_position),
              manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, textsize = 5) +
  scale_color_manual(values = condition_surface_colors, guide = guide_legend(nrow = 3, byrow = FALSE)) +
  scale_shape_manual(values = surface_shapes, breaks = surface_order) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  guides(shape = "none") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), panel.grid.minor = element_blank()) +
  labs(y = "QCLog_CmbCo intensity", x = "")

sessionInfo()