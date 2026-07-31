## =============================================================================
## UKA (Upstream Kinase Analysis) — Kinase Activity Heatmaps
## =============================================================================

# ---- Packages ----------------------------------------------------------------
library(dplyr)
library(readr)
library(pheatmap)

# ---- Paths --------------------------------------------------------------------
# Edit `data_dir` to point at the folder containing the four UKA export CSVs.
# Using file.path() keeps this portable across operating systems.
data_dir <- "YOUR/DATA/REPOSITORY/"

path_hard_vs_soft    <- file.path(data_dir, "1_UKA_Hard vs Soft Run1-2.csv")
path_stretch_vs_soft <- file.path(data_dir, "3_UKA_Stretch vs Soft Run3-4.csv")
path_hardsoft_geno   <- file.path(data_dir, "2_UKA_LMNA vs WT run1-2.csv")
path_stretch_geno    <- file.path(data_dir, "4_UKA_LMNA vs WT run3-4.csv")

score_threshold_substrate <- 1.3
score_threshold_genotype  <- 1.3

# ---- Shared helper functions ---------------------------------------------------
extract_contrast <- function(df, contrast_col, pattern, id_cols, suffix) {
  df %>%
    filter(grepl(pattern, .data[[contrast_col]])) %>%
    rename_with(~ paste0(., suffix), .cols = -all_of(id_cols))
}
zero_below_threshold <- function(df, stat_col, score_col, threshold) {
  df[[stat_col]] <- ifelse(df[[score_col]] < threshold, 0, df[[stat_col]])
  df
}
plot_kinase_heatmap <- function(mat, col_order = NULL, flip_cols = NULL,
                                 gaps_row = NULL, main = "") {
  if (!is.null(col_order)) mat <- mat[, col_order, drop = FALSE]
  if (!is.null(flip_cols)) {
    flip_idx <- grepl(flip_cols, colnames(mat))
    mat[, flip_idx] <- -mat[, flip_idx]
  }

  palette <- c(
    colorRampPalette(c("blue", "white"))(25),
    "grey",
    colorRampPalette(c("white", "red"))(25)
  )
  breaks <- c(seq(-1, -0.001, length.out = 26), seq(0, 15, length.out = 26)[-1])

  pheatmap(
    mat               = t(mat),
    cluster_rows      = FALSE,
    cluster_cols      = TRUE,
    display_numbers   = FALSE,
    color             = palette,
    breaks            = breaks,
    labels_col        = rownames(mat),
    fontsize_col      = 10,
    fontsize_row      = 8,
    angle_col         = 90,
    gaps_row          = gaps_row,
    main              = main
  )
}

# =============================================================================
## 1. Substrate mechanics comparison heatmap
# =============================================================================

hard_vs_soft    <- read_csv(path_hard_vs_soft, show_col_types = FALSE)
stretch_vs_soft <- read_csv(path_stretch_vs_soft, show_col_types = FALSE)

id_cols_hs <- c("UKA_app.UKA.Kinase Name", "UKA_app.UKA.Kinase Uniprot ID")
id_cols_ss <- c("UKA_app1.UKA1.Kinase Name", "UKA_app1.UKA1.Kinase Uniprot ID")

wt_hard_vs_soft   <- extract_contrast(hard_vs_soft, "UKA_app.UKA.Sgroup_contrast",
                                       "C_WT_Hard vs C_Soft", id_cols_hs, "_WT_HARDvsSOFT")
lmna_hard_vs_soft <- extract_contrast(hard_vs_soft, "UKA_app.UKA.Sgroup_contrast",
                                       "LMNA_Hard vs C_Soft", id_cols_hs, "_LMNA_HARDvsSOFT")
wt_stretch_vs_soft   <- extract_contrast(stretch_vs_soft, "UKA_app1.UKA1.Sgroup_contrast",
                                          "C_WT_Stretch vs C_Soft", id_cols_ss, "_WT_STRETCHvsSOFT")
lmna_stretch_vs_soft <- extract_contrast(stretch_vs_soft, "UKA_app1.UKA1.Sgroup_contrast",
                                          "LMNA_Stretch vs C_Soft", id_cols_ss, "_LMNA_STRETCHvsSOFT")

hardsoft    <- merge(wt_hard_vs_soft, lmna_hard_vs_soft, by = "UKA_app.UKA.Kinase Name")
stretchsoft <- merge(wt_stretch_vs_soft, lmna_stretch_vs_soft, by = "UKA_app1.UKA1.Kinase Name") %>%
  rename(`UKA_app.UKA.Kinase Name` = `UKA_app1.UKA1.Kinase Name`)

substrate_uka <- merge(hardsoft, stretchsoft, by = "UKA_app.UKA.Kinase Name") %>%
  zero_below_threshold("UKA_app1.UKA1.Mean Kinase Statistic_LMNA_STRETCHvsSOFT",
                        "UKA_app1.UKA1.Median Final score_LMNA_STRETCHvsSOFT", score_threshold_substrate) %>%
  zero_below_threshold("UKA_app1.UKA1.Mean Kinase Statistic_WT_STRETCHvsSOFT",
                        "UKA_app1.UKA1.Median Final score_WT_STRETCHvsSOFT", score_threshold_substrate) %>%
  zero_below_threshold("UKA_app.UKA.Mean Kinase Statistic_LMNA_HARDvsSOFT",
                        "UKA_app.UKA.Median Final score_LMNA_HARDvsSOFT", score_threshold_substrate) %>%
  zero_below_threshold("UKA_app.UKA.Mean Kinase Statistic_WT_HARDvsSOFT",
                        "UKA_app.UKA.Median Final score_WT_HARDvsSOFT", score_threshold_substrate)

heatmap_substrate <- substrate_uka %>%
  select(`UKA_app.UKA.Kinase Name`, contains("Mean Kinase Statistic")) %>%
  rename(
    `hfLMNA_R377L \nSoft strained vs Soft unstrained`    = `UKA_app1.UKA1.Mean Kinase Statistic_LMNA_STRETCHvsSOFT`,
    `hfLMNA_WT \nSoft strained vs Soft unstrained`       = `UKA_app1.UKA1.Mean Kinase Statistic_WT_STRETCHvsSOFT`,
    `hfLMNA_R377L \nRigid unstrained vs Soft unstrained` = `UKA_app.UKA.Mean Kinase Statistic_LMNA_HARDvsSOFT`,
    `hfLMNA_WT \nRigid unstrained vs Soft unstrained`    = `UKA_app.UKA.Mean Kinase Statistic_WT_HARDvsSOFT`
  ) %>%
  as.data.frame()

rownames(heatmap_substrate) <- make.names(heatmap_substrate$`UKA_app.UKA.Kinase Name`)
heatmap_substrate <- heatmap_substrate %>%
  select(-`UKA_app.UKA.Kinase Name`) %>%
  as.matrix()

substrate_col_order <- c(
  "hfLMNA_WT \nRigid unstrained vs Soft unstrained",
  "hfLMNA_WT \nSoft strained vs Soft unstrained",
  "hfLMNA_R377L \nRigid unstrained vs Soft unstrained",
  "hfLMNA_R377L \nSoft strained vs Soft unstrained"
)

plot_kinase_heatmap(
  heatmap_substrate,
  col_order = substrate_col_order,
  flip_cols = "Rigid",
  gaps_row  = 2,
  main      = "Kinase activity"
)

# =============================================================================
## 2. Genotype comparison heatmap
# =============================================================================

hardsoft_geno <- read_csv(path_hardsoft_geno, show_col_types = FALSE)
stretch_geno  <- read_csv(path_stretch_geno, show_col_types = FALSE)

id_cols_hsg <- c("UKA_app0.UKA0.Kinase Name", "UKA_app0.UKA0.Kinase Uniprot ID")
id_cols_stg <- c("UKA_app10.UKA10.Kinase Name", "UKA_app10.UKA10.Kinase Uniprot ID")

hard_geno    <- extract_contrast(hardsoft_geno, "UKA_app0.UKA0.Sgroup_contrast",
                                  "Hard_LMNA vs C_WT", id_cols_hsg, "_Hard_LMNAvsWT")
soft_geno    <- extract_contrast(hardsoft_geno, "UKA_app0.UKA0.Sgroup_contrast",
                                  "C_Soft_LMNA vs C_WT", id_cols_hsg, "_Soft_LMNAvsWT")
stretch_geno <- extract_contrast(stretch_geno, "UKA_app10.UKA10.Sgroup_contrast",
                                  "Stretch_LMNA vs C_WT", id_cols_stg, "_Stretch_LMNAvsWT")

hard_soft_geno <- merge(hard_geno, soft_geno, by = "UKA_app0.UKA0.Kinase Name")
stretch_geno   <- stretch_geno %>%
  rename(`UKA_app0.UKA0.Kinase Name` = `UKA_app10.UKA10.Kinase Name`)

genotype_uka <- merge(hard_soft_geno, stretch_geno, by = "UKA_app0.UKA0.Kinase Name") %>%
  zero_below_threshold("UKA_app0.UKA0.Mean Kinase Statistic_Hard_LMNAvsWT",
                        "UKA_app0.UKA0.Median Final score_Hard_LMNAvsWT", score_threshold_genotype) %>%
  zero_below_threshold("UKA_app0.UKA0.Mean Kinase Statistic_Soft_LMNAvsWT",
                        "UKA_app0.UKA0.Median Final score_Soft_LMNAvsWT", score_threshold_genotype) %>%
  zero_below_threshold("UKA_app10.UKA10.Mean Kinase Statistic_Stretch_LMNAvsWT",
                        "UKA_app10.UKA10.Median Final score_Stretch_LMNAvsWT", score_threshold_genotype)

heatmap_genotype <- genotype_uka %>%
  select(`UKA_app0.UKA0.Kinase Name`, contains("Mean Kinase Statistic")) %>%
  rename(
    `Soft unstrained`  = `UKA_app0.UKA0.Mean Kinase Statistic_Soft_LMNAvsWT`,
    `Rigid unstrained` = `UKA_app0.UKA0.Mean Kinase Statistic_Hard_LMNAvsWT`,
    `Soft strained`    = `UKA_app10.UKA10.Mean Kinase Statistic_Stretch_LMNAvsWT`
  ) %>%
  as.data.frame()

rownames(heatmap_genotype) <- make.names(heatmap_genotype$`UKA_app0.UKA0.Kinase Name`)
heatmap_genotype <- heatmap_genotype %>%
  select(-`UKA_app0.UKA0.Kinase Name`) %>%
  as.matrix()

plot_kinase_heatmap(
  heatmap_genotype,
  main = "Kinase activity (hfLMNA_R377L vs hfLMNA_WT)"
)
