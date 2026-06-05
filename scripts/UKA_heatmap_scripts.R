### UKA analysis ###

# Packages
library(dplyr)
library(pheatmap)
library(readxl)
library(ggplot2)
library(readr)

# Load data ### ADJUST PATH ###
Hard_vs_soft <- read_csv("YOUR/LOCAL/PATH/1_UKA_Hard vs Soft Run1-2.csv")
Stretch_vs_soft <- read_csv("YOUR/LOCAL/PATH/3_UKA_Stretch vs Soft Run3-4.csv")

WT_hard_vs_soft <- Hard_vs_soft[grepl("C_WT_Hard vs C_Soft", Hard_vs_soft$UKA_app.UKA.Sgroup_contrast),]
LMNA_hard_vs_soft <- Hard_vs_soft[grepl("LMNA_Hard vs C_Soft", Hard_vs_soft$UKA_app.UKA.Sgroup_contrast),]
WT_stretch_vs_soft <- Stretch_vs_soft[grepl("C_WT_Stretch vs C_Soft", Stretch_vs_soft$UKA_app1.UKA1.Sgroup_contrast),]
LMNA_stretch_vs_soft <- Stretch_vs_soft[grepl("LMNA_Stretch vs C_Soft", Stretch_vs_soft$UKA_app1.UKA1.Sgroup_contrast),]

WT_hard_vs_soft <- WT_hard_vs_soft %>%
  rename_with(~ paste0(., "_WT_HARDvsSOFT"), -c(`UKA_app.UKA.Kinase Name`, `UKA_app.UKA.Kinase Uniprot ID`))
LMNA_hard_vs_soft <- LMNA_hard_vs_soft %>%
  rename_with(~ paste0(., "_LMNA_HARDvsSOFT"), -c(`UKA_app.UKA.Kinase Name`, `UKA_app.UKA.Kinase Uniprot ID`))
WT_stretch_vs_soft <- WT_stretch_vs_soft %>%
  rename_with(~ paste0(., "_WT_STRETCHvsSOFT"), -c(`UKA_app1.UKA1.Kinase Name`, `UKA_app1.UKA1.Kinase Uniprot ID`))
LMNA_stretch_vs_soft <- LMNA_stretch_vs_soft %>%
  rename_with(~ paste0(., "_LMNA_STRETCHvsSOFT"), -c(`UKA_app1.UKA1.Kinase Name`, `UKA_app1.UKA1.Kinase Uniprot ID`))

hardsoft <- merge(WT_hard_vs_soft, LMNA_hard_vs_soft, by = "UKA_app.UKA.Kinase Name")
stretchsoft <- merge(WT_stretch_vs_soft, LMNA_stretch_vs_soft, by = "UKA_app1.UKA1.Kinase Name")
names(stretchsoft)[names(stretchsoft) == "UKA_app1.UKA1.Kinase Name"] <- "UKA_app.UKA.Kinase Name"
sub_UKA <- merge(hardsoft, stretchsoft, by = "UKA_app.UKA.Kinase Name")
sub_UKA <- sub_UKA %>%
  mutate(
    `UKA_app1.UKA1.Mean Kinase Statistic_LMNA_STRETCHvsSOFT` = 
      ifelse(`UKA_app1.UKA1.Median Final score_LMNA_STRETCHvsSOFT` < 1.3, 0, `UKA_app1.UKA1.Mean Kinase Statistic_LMNA_STRETCHvsSOFT`),
    
    `UKA_app1.UKA1.Mean Kinase Statistic_WT_STRETCHvsSOFT` = 
      ifelse(`UKA_app1.UKA1.Median Final score_WT_STRETCHvsSOFT` < 1.3, 0, `UKA_app1.UKA1.Mean Kinase Statistic_WT_STRETCHvsSOFT`),
    
    
    `UKA_app.UKA.Mean Kinase Statistic_LMNA_HARDvsSOFT` = 
      ifelse(`UKA_app.UKA.Median Final score_LMNA_HARDvsSOFT` < 1.3, 0, `UKA_app.UKA.Mean Kinase Statistic_LMNA_HARDvsSOFT`),
    
    `UKA_app.UKA.Mean Kinase Statistic_WT_HARDvsSOFT` = 
      ifelse(`UKA_app.UKA.Median Final score_WT_HARDvsSOFT` < 1.3, 0, `UKA_app.UKA.Mean Kinase Statistic_WT_HARDvsSOFT`)
  )

heatmap_SUB <- sub_UKA %>%
  dplyr::select(
    `UKA_app.UKA.Kinase Name`,
    contains("Mean Kinase Statistic")
  ) %>%
  dplyr::rename(
    `hfLMNA_R377L \nSoft strained vs Soft unstrained`    = `UKA_app1.UKA1.Mean Kinase Statistic_LMNA_STRETCHvsSOFT`,
    `hfLMNA_WT \nSoft strained vs Soft unstrained` = `UKA_app1.UKA1.Mean Kinase Statistic_WT_STRETCHvsSOFT`,
    `hfLMNA_R377L \nRigid unstrained vs Soft unstrained` = `UKA_app.UKA.Mean Kinase Statistic_LMNA_HARDvsSOFT`,
    `hfLMNA_WT \nRigid unstrained vs Soft unstrained` = `UKA_app.UKA.Mean Kinase Statistic_WT_HARDvsSOFT`
  ) %>%
  as.data.frame()
rownames(heatmap_SUB) <- make.names(heatmap_SUB$`UKA_app.UKA.Kinase Name`)
heatmap_SUB <- heatmap_SUB %>%
  dplyr::select(-`UKA_app.UKA.Kinase Name`)

heatmap_SUB <- as.matrix(heatmap_SUB)

# UKA substrate comparison Heatmap
my_palette <- c(
  colorRampPalette(c("blue", "white"))(25),
  "grey",  
  colorRampPalette(c("white", "red"))(25)
)
breaks <- c(seq(-1, -0.001, length.out = 26), seq(0, 15, length.out = 26)[-1])
desired_order <- c(
  "hfLMNA_WT \nRigid unstrained vs Soft unstrained",
  "hfLMNA_WT \nSoft strained vs Soft unstrained",
  "hfLMNA_R377L \nRigid unstrained vs Soft unstrained",
  "hfLMNA_R377L \nSoft strained vs Soft unstrained"
)
heatmap_SUB <- heatmap_SUB[, desired_order]
heatmap_SUB[, grepl("Rigid", colnames(heatmap_SUB))] <-
  -heatmap_SUB[, grepl("Rigid", colnames(heatmap_SUB))]

pheatmap(
  mat = t(heatmap_SUB),
  cluster_rows = FALSE, 
  cluster_cols = TRUE, 
  display_numbers = FALSE,
  color = my_palette,
  breaks = breaks,
  labels_col = rownames(heatmap_SUB),
  fontsize_col = 10,
  fontsize_row = 8,
  angle_col = 90,
  gaps_row = 2,
  main = "Kinase activity"
)

# UKA Genotype comparison Heatmap
HardSoft_LMNAvsWT <- read_csv("YOUR/LOCAL/PATH/2_UKA_ LMNA vs WT run1-2.csv")
StretchSoft_LMNAvsWT <- read_csv("YOUR/LOCAL/PATH/4_UKA_LMNA vs WT run3-4.csv")

Hard_LMNAvsWT <- HardSoft_LMNAvsWT[grepl("Hard_LMNA vs C_WT", HardSoft_LMNAvsWT$UKA_app0.UKA0.Sgroup_contrast),]
Soft_LMNAvsWT <- HardSoft_LMNAvsWT[grepl("C_Soft_LMNA vs C_WT", HardSoft_LMNAvsWT$UKA_app0.UKA0.Sgroup_contrast),]
Stretch_LMNAvsWT <- StretchSoft_LMNAvsWT[grepl("Stretch_LMNA vs C_WT", StretchSoft_LMNAvsWT$UKA_app10.UKA10.Sgroup_contrast),]

Hard_LMNAvsWT <- Hard_LMNAvsWT %>%
  rename_with(~ paste0(., "_Hard_LMNAvsWT"), -c(`UKA_app0.UKA0.Kinase Name`, `UKA_app0.UKA0.Kinase Uniprot ID`))
Soft_LMNAvsWT <- Soft_LMNAvsWT %>%
  rename_with(~ paste0(., "_Soft_LMNAvsWT"), -c(`UKA_app0.UKA0.Kinase Name`, `UKA_app0.UKA0.Kinase Uniprot ID`))
Stretch_LMNAvsWT <- Stretch_LMNAvsWT %>%
  rename_with(~ paste0(., "_Stretch_LMNAvsWT"), -c(`UKA_app10.UKA10.Kinase Name`, `UKA_app10.UKA10.Kinase Uniprot ID`))

HardSoft <- merge(Hard_LMNAvsWT, Soft_LMNAvsWT, by = "UKA_app0.UKA0.Kinase Name")
names(Stretch_LMNAvsWT)[names(Stretch_LMNAvsWT) == "UKA_app10.UKA10.Kinase Name"] <- "UKA_app0.UKA0.Kinase Name"

sub_UKA <- merge(HardSoft, Stretch_LMNAvsWT, by = "UKA_app0.UKA0.Kinase Name")
sub_UKA <- sub_UKA %>%
  mutate(
    `UKA_app0.UKA0.Mean Kinase Statistic_Hard_LMNAvsWT` = 
      ifelse(`UKA_app0.UKA0.Median Final score_Hard_LMNAvsWT` < 1.2, 0, `UKA_app0.UKA0.Mean Kinase Statistic_Hard_LMNAvsWT`),
    
    `UKA_app0.UKA0.Mean Kinase Statistic_Soft_LMNAvsWT` = 
      ifelse(`UKA_app0.UKA0.Median Final score_Soft_LMNAvsWT` < 1.2, 0, `UKA_app0.UKA0.Mean Kinase Statistic_Soft_LMNAvsWT`),
    
    `UKA_app10.UKA10.Mean Kinase Statistic_Stretch_LMNAvsWT` = 
      ifelse(`UKA_app10.UKA10.Median Final score_Stretch_LMNAvsWT` < 1.2, 0, `UKA_app10.UKA10.Mean Kinase Statistic_Stretch_LMNAvsWT`)
  )

heatmap_SUB <- sub_UKA %>%
  dplyr::select(
    `UKA_app0.UKA0.Kinase Name`,
    contains("Mean Kinase Statistic")
  ) %>%
  dplyr::rename(
    `Soft unstrained`    = `UKA_app0.UKA0.Mean Kinase Statistic_Soft_LMNAvsWT`,
    `Rigid unstrained` = `UKA_app0.UKA0.Mean Kinase Statistic_Hard_LMNAvsWT`,
    `Soft strained` = `UKA_app10.UKA10.Mean Kinase Statistic_Stretch_LMNAvsWT`
  ) %>%
  as.data.frame()
rownames(heatmap_SUB) <- make.names(heatmap_SUB$`UKA_app0.UKA0.Kinase Name`)
heatmap_SUB <- heatmap_SUB %>%
  dplyr::select(-`UKA_app0.UKA0.Kinase Name`)
heatmap_SUB <- as.matrix(heatmap_SUB)

# UKA genotype comparison heatmap
my_palette <- c(
  colorRampPalette(c("blue", "white"))(25),
  "grey",  
  colorRampPalette(c("white", "red"))(25)
)
breaks <- c(seq(-1, -0.001, length.out = 26), seq(0, 15, length.out = 26)[-1])
pheatmap(
  mat = t(heatmap_SUB),
  cluster_rows = FALSE, 
  cluster_cols = TRUE, 
  display_numbers = FALSE,
  color = my_palette,
  breaks = breaks,
  labels_col = rownames(heatmap_SUB),
  fontsize_col = 10,
  fontsize_row = 8,
  angle_col = 90,
  main = "Kinase activity (hfLMNA_R377L vs hfLMNA_WT)"
)
