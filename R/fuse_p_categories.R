library(tidyverse)
library(ggcorrplot)
library(plotly)
library(kableExtra)
library(export)
library(RColorBrewer)
library(ggpubr)

rm(list = ls())

data_path <- "outputs/data"
figure_path <- "outputs/figures/"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")
source("R/util_lib.R")

gene_target <- "BRCA1"
cancer_filenames <- getCancerFiles("data/BRCA1", gene_target)
filename <- "20231201_BRCA1_Fuse_PVal.pptx"
file_path <- paste0(figure_path, filename)
size_col <- "beta"

mave_data <- readMaveData("data/mave_data_brn_v3.csv")

# Breast Cancer ----
data <- get_genebass_mave_data_with_pval(
  gene_target,
  mave_data,
  "BreastCancer",
  cancer_filenames
)

data2 <- data %>%
  filter((pval >= 0.1) & (FUSE_score < -1.5)) %>%
  select(variant_id, ClinVarLabelP, FUSE_score, 
         lower_bound, upper_bound, is_outlier, is_lower_outlier, is_upper_outlier,
         is_outlier_label, row_id)

title_txt <- paste(gene_target, "Breast Cancer",
  "FUSE Score",
  "Variants by P Value Category",
  sep = ","
)
xlabel <- "P Value Category"
ylabel <- "FUSE Score"


pt <- get_violin_box_FUSE_score_by_pval_category(data,
  title_txt = title_txt,
  size_col = size_col,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 5,
  ybreaks =  seq(-3, 3, by = 0.5)
)
pt <- pt + stat_compare_means(label.y = 2.5)
pt

# pt <- cowplot::plot_grid(pt1, pt2, nrow = 2)
graph2ppt(pt, file_path, width = 7, height = 7, append = TRUE)


# Colorectal Cancer ----
data <- get_genebass_mave_data_with_pval(
  gene_target,
  mave_data,
  "ColorectalCancer",
  cancer_filenames
)
title_txt <- paste(gene_target, "Colorectal Cancer",
  "FUSE Score",
  "Variants by P Value Category",
  sep = ","
)
xlabel <- "P Value Category"
ylabel <- "FUSE Score"

pt <- get_violin_box_FUSE_score_by_pval_category(data,
  title_txt = title_txt,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 2
)
pt <- pt + stat_compare_means(label.y = 3)
pt
graph2ppt(pt, file_path, width = 7, height = 7, append = TRUE)


# Lung Cancer ----
data <- get_genebass_mave_data_with_pval(
  gene_target,
  mave_data,
  "LungCancer",
  cancer_filenames
) %>%
  filter(!as.character(pval_category) == "[0, 0.01]")


title_txt <- paste(gene_target, "Lung Cancer",
  "FUSE Score",
  "Variants by P Value Category",
  sep = ","
)
xlabel <- "P Value Category"
ylabel <- "FUSE Score"

pt <- get_violin_box_FUSE_score_by_pval_category(data,
  title_txt = title_txt,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 2
)
pt <- pt + stat_compare_means(label.y = 3)
pt
graph2ppt(pt, file_path, width = 7, height = 7, append = TRUE)

# Prostate Cancer ----

data <- get_genebass_mave_data_with_pval(
  gene_target,
  mave_data,
  "ProstateCancer",
  cancer_filenames
)

title_txt <- paste(gene_target, "Prostate Cancer",
  "FUSE Score",
  "Variants by P Value Category",
  sep = ","
)
xlabel <- "P Value Category"
ylabel <- "FUSE Score"

pt <- get_violin_box_FUSE_score_by_pval_category(data,
  title_txt = title_txt,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 2
)
pt <- pt + stat_compare_means(label.y = -3)
pt
graph2ppt(pt, file_path, width = 7, height = 7, append = TRUE)


# Rare Cancer ----

data <- get_genebass_mave_data_with_pval(
  gene_target,
  mave_data,
  "ProstateCancer",
  cancer_filenames
)

title_txt <- paste(gene_target, "Prostate Cancer",
                   "FUSE Score",
                   "Variants by P Value Category",
                   sep = ","
)
xlabel <- "P Value Category"
ylabel <- "FUSE Score"

pt <- get_violin_box_FUSE_score_by_pval_category(data,
                                                 title_txt = title_txt,
                                                 xlabel = xlabel,
                                                 ylabel = ylabel,
                                                 annotate_flag = TRUE,
                                                 annotate_text_size = 2
)
pt <- pt + stat_compare_means(label.y = -3)
pt
graph2ppt(pt, file_path, width = 7, height = 7, append = TRUE)
