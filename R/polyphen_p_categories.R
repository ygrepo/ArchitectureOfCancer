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

#Get the gene name from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No gene name provided. Usage: Rscript polyphen_p_categories.R <gene_name>", call. = FALSE)
}
gene <- args[1]

ppt_filename <- paste0("20240121_", gene, "_Polyphen_PVal.pptx")
size_col <- "BETA"
var_label_col <- "VarLabel"
xlabel <- "P Value Category"
ylabel <- "Polyphen Score"
percent_IQRValue <- 0.25

df <- read_csv_polyphen_df(gene)


for (cancer_type in unique(df$description)) {
  if (is.na(cancer_type) | (cancer_type == "")) {
    next
  }
  print(cancer_type)
  cancer_df <- process_polyphen_data_with_pval(df,
    description_val = cancer_type,
    percent_IQR = percent_IQRValue
  )


  title_txt <- paste(gene, cancer_type,
    "Polyphen Score",
    "Variants by P Value Category",
    sep = ","
  )

  pt <- get_violin_box_polyphen_score_by_pval_category(df = cancer_df,
    title_txt = title_txt,
    size_col = size_col,
    var_label_col = var_label_col,
    xlabel = xlabel,
    ylabel = ylabel,
    title_font_size = 10,
    annotate_flag = TRUE,
    annotate_text_size = 5,
    ybreaks = seq(-1, 1, by = 0.25)
  )
  pt

  save_to_ppt(filename = ppt_filename, pt)
}

