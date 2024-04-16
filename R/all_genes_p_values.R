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



genes <- read_csv_input_data(filename = "gene_list.txt", colnames = c("Gene"))
result_df <- tibble()
for (gene in genes$Gene) {
  print(paste0("Gene:", gene))
  df <- read_csv_polyphen_df(gene)
  print(unique(df$description))
  result_df <- bind_rows(result_df, df)
}


result_df <- result_df %>%
  dplyr::filter(!is.na(description)) %>%
  dplyr::filter((description != "") & (description != "NA"))

cancer_types <- unique(result_df$description)
cancer_types
genes <- unique(result_df$gene)
genes
genes <- c("MLH1", "MSH2", "MSH6","PMS2")
#genes <- c("BRCA1", "BRCA2", "PALB2")
significan_test_df <- tibble()
for (cancer_type in cancer_types) {
  for (gene in genes) {
    print(paste0("Cancer:", cancer_type, ",Gene:", gene))
    df <- result_df %>%
      create_pval_category(
        gene_val = gene,
        description_val = cancer_type
      )
    if (length(df$pval_category) == 0) {
      next
    }
    df <- df %>%
      perform_polyphen_score_signif_assoc_by_p_val_category(
        cancer_type = cancer_type,
        gene_val = gene
      )

    significan_test_df <- bind_rows(significan_test_df, df)
  }
}

significan_test_df <- significan_test_df %>%
  dplyr::filter(!is.na(p_value))


cancer_types <- unique(significan_test_df$cancer_type)
cancer_types
cancer_types <- c("Prostate cancer")
#cancer_types <- c("Breast cancer")

significan_test_adj_pval_df <- tibble()
for (cancer_type in cancer_types) {
  print(paste0("Cancer type:", cancer_type))
  df <- significan_test_df %>%
    dplyr::filter(cancer_type == cancer_type) %>%
    dplyr::mutate(pval_cat = paste0(group1, group2)) %>%
    dplyr::group_by(pval_cat) %>%
    dplyr::arrange(pval_cat) %>%
    dplyr::mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
    dplyr::mutate(LOG10_adj_pval = -log10(adjusted_p_value))
  significan_test_adj_pval_df <- bind_rows(significan_test_adj_pval_df, df)
}

significan_test_adj_pval_df <- significan_test_adj_pval_df %>%
  dplyr::distinct()

unique(significan_test_adj_pval_df$pval_cat)
df2 <- significan_test_adj_pval_df %>%
  dplyr::filter(cancer_type == "Prostate cancer") 
#  dplyr::filter(cancer_type == "Breast cancer") 
# %>%
#   dplyr::filter(adjusted_p_value < 0.5)
unique(df2$pval_cat)
df2 %>%
  filter(gene == "MLH1") %>%
  filter(pval_cat == "[0, 0.01]]0.01, 0.1]")

filename <- "20240415_Polyphen_Genes_Significance_Tests_ColorectalCancer_type.csv"
#filename <- "20240302_Polyphen_Genes_Significance_Tests_BreastCancer_type.csv"
write_csv_data(df = significan_test_adj_pval_df, filename = filename)


title_txt <- "Colorectal Cancer, Polyphen Score Significance Association Test"
#title_txt <- "Breast Cancer, Polyphen Score Significance Association Test"
pt <- get_pval_category_adj_pval_scatter_plot(df = df2,
                                        title_txt = title_txt,
                                        max_overlaps_val = 50
                                        )
pt
ppt_filename <- "20240415_Polyphen_Significance_Association_Test.pptx"
#ppt_filename <- "20240302_Polyphen_Significance_Association_Test.pptx"
save_to_ppt(filename = ppt_filename, pt)

title_txt <- "All Cancers, Polyphen Score Significance Association Test"
pt <- get_pval_category_adj_pval_scatter_plot(df = significan_test_adj_pval_df,
                                              title_txt = title_txt,
                                              max_overlaps_val = 50
)
pt
save_to_ppt(filename = ppt_filename, pt)

