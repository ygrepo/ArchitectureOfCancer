library(tidyverse)
library(dplyr)
library(ggcorrplot)
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

gene_target <- "BRCA1"
ppt_filename <- "20240117_BRCA1_Polyphen_PVal.pptx"
size_col <- "beta"
xlabel <- "P Value Category"
ylabel <- "Polyphen Score"

filename <- "outputs/data/BRCA1_breast_cancer_polyphen.csv"
data <- as_tibble(read.table(filename,
  header = TRUE, sep = " "
))


# Breast Cancer ----
cancer_type <- "Breast cancer"
title_cancer_type <- "Breast Cancer"
df <- get_genebass_polyphen_data_with_pval(data,
  description_val = cancer_type,
  percent_IQR = 0.25
)

title_txt <- paste(gene_target, title_cancer_type,
  "Polyphen Score",
  "Variants by P Value Category",
  sep = ","
)

pt <- get_violin_box_polyphen_score_by_pval_category(df,
  title_txt = title_txt,
  size_col = size_col,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 5,
  ybreaks = seq(-1, 1, by = 0.25)
)
pt

save_to_ppt(filename = ppt_filename, pt)

# Lung cancer ----
cancer_type <- "Lung cancer"
title_cancer_type <- "Lung Cancer"
df <- get_genebass_polyphen_data_with_pval(data,
  description_val = cancer_type,
  percent_IQR = 0.25
)

title_txt <- paste(gene_target, title_cancer_type,
  "Polyphen Score",
  "Variants by P Value Category",
  sep = ","
)

pt <- get_violin_box_polyphen_score_by_pval_category(df,
  title_txt = title_txt,
  size_col = size_col,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 5,
  ybreaks = seq(-1, 1, by = 0.25)
)
pt

save_to_ppt(filename = ppt_filename, pt)

# Prostate cancer ----
cancer_type <- "Prostate cancer"
title_cancer_type <- "Prostate Cancer"
df <- get_genebass_polyphen_data_with_pval(data,
  description_val = cancer_type,
  percent_IQR = 0.25
)

title_txt <- paste(gene_target, title_cancer_type,
  "Polyphen Score",
  "Variants by P Value Category",
  sep = ","
)

pt <- get_violin_box_polyphen_score_by_pval_category(df,
  title_txt = title_txt,
  size_col = size_col,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 5,
  ybreaks = seq(-1, 1, by = 0.25)
)
pt

save_to_ppt(filename = ppt_filename, pt)

# Bowel cancer ----
cancer_type <- "Bowel cancer in the colon or rectum"
title_cancer_type <- "Bowel Cancer"
df <- get_genebass_polyphen_data_with_pval(data,
  description_val = cancer_type,
  percent_IQR = 0.25
)

title_txt <- paste(gene_target, title_cancer_type,
  "Polyphen Score",
  "Variants by P Value Category",
  sep = ","
)

pt <- get_violin_box_polyphen_score_by_pval_category(df,
  title_txt = title_txt,
  size_col = size_col,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 5,
  ybreaks = seq(-1, 1, by = 0.25)
)
pt

save_to_ppt(filename = ppt_filename, pt)

# Rare cancer ----
cancer_type <- "Rare cancer"
title_cancer_type <- "Rare Cancer"
df <- get_genebass_polyphen_data_with_pval(data,
  description_val = cancer_type,
  percent_IQR = 0.25
)

title_txt <- paste(gene_target, title_cancer_type,
  "Polyphen Score",
  "Variants by P Value Category",
  sep = ","
)

pt <- get_violin_box_polyphen_score_by_pval_category(df,
  title_txt = title_txt,
  size_col = size_col,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 5,
  ybreaks = seq(-1, 1, by = 0.25)
)
pt

save_to_ppt(filename = ppt_filename, pt)

# "Cancer code, self-reported ----
cancer_type <- "Cancer code, self-reported"
title_cancer_type <- "Cancer code, self-reported"
df <- get_genebass_polyphen_data_with_pval(data,
  description_val = cancer_type,
  percent_IQR = 0.25
)

title_txt <- paste(gene_target, title_cancer_type,
  "Polyphen Score",
  "Variants by P Value Category",
  sep = ","
)

pt <- get_violin_box_polyphen_score_by_pval_category(df,
  title_txt = title_txt,
  size_col = size_col,
  xlabel = xlabel,
  ylabel = ylabel,
  annotate_flag = TRUE,
  annotate_text_size = 5,
  ybreaks = seq(-1, 1, by = 0.25)
)
pt

save_to_ppt(filename = ppt_filename, pt)
