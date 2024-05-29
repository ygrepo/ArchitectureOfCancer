library(tidyverse)
library(dplyr)
library(ggcorrplot)
library(ggrepel)
library(plotly)
library(kableExtra)
library(export)
library(RColorBrewer)
library(ggpubr)

rm(list = ls())


setwd("~/github/ArchitectureOfCancer/")

source("R/io_utils.R")
source("R/util_lib.R")
source("R/plot_lib.R")

tissue <- "BreastCancer"
size_col <- "beta"
xlabel <- "P Value Category"
ylabel <- "Polyphen Score"


# Breast Cancer ----
cancer_type <- "Breast cancer"
gene_val <- "PALB2"
tissue <- "BreastCancer"

cancer_filenames <- getCancerFiles(gene_val)
genebass_df <- read_genebass_data_with_variant_process(
  cancer_filenames = cancer_filenames,
  tissue = tissue
)

signifPVal <- (-log10(0.05))
df_cancer <- read_csv_polyphen_df(gene_val) %>%
  dplyr::filter((description == cancer_type) & (gene == gene_val)) %>%
  na.omit((polyphen_score)) %>%
  na.omit(Pvalue) %>%
  dplyr::mutate(Log10PVal = -log10(Pvalue)) %>%
  dplyr::mutate(correct_AF = ifelse(AF >= 0.5, 1 - AF, AF)) %>%
  dplyr::mutate(LOGAF = -log10(correct_AF))

df_cancer <- df_cancer %>%
  dplyr::left_join(genebass_df) %>%
  dplyr::mutate(variantIdSign = case_when(
    Log10PVal <= signifPVal ~ "",
    TRUE ~ shgvsp
  ))

summary(df_cancer$LOGAF)
title_txt <- "Breast Cancer, PALB2, Polyphen Score vs. -LOG10(PValue)"
alpha <- 1
pt <- get_polyphen_score_genebass_pvalue(
  df = df_cancer,
  title_txt,
  alpha = alpha,
  title_font_size = 24,
  x_y_font_size = 24,
  legend_title_font_size = 14,
  legend_text_font_size = 8,
  stats_font_size = 6
)
pt

filename <- "20240519_BreastCancer_PALB2_Polyphen_LOG10PVal.jpg"
save_plot(
  filename = filename,
  pt = pt
)
