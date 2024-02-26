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


# cancer_type <- "Breast cancer"

# df1 <- read_csv_polyphen_df("APC") %>%
#   dplyr::filter(description == cancer_type)
# df1
#
# df2 <- read_csv_polyphen_df("BRCA1") %>%
#   dplyr::filter(description == cancer_type)
# unique(df1$polyphen_score) %>% sort()
# unique(df2$polyphen_score) %>% sort()
#
# unique(df1$Pvalue)
#
# result_df <- bind_rows(df1, df2)
#
# df11 <- result_df %>%
#   create_pval_category(
#     gene_val = "APC",
#     description_val = cancer_type
#   )
# unique(df11$pval_category)
# df11

# t1 <- perform_polyphen_score_signif_assoc_by_p_val_category(
#   df = df11,
#   cancer_type = cancer_type,
#   gene_val = "APC"
# )
# t1
#
# df21 <- result_df %>%
#   create_pval_category(
#     gene_val = "BRCA1",
#     description_val = cancer_type
#   )
#
# unique(df21$pval_category)

# df3 <- read_csv_polyphen_df("TSC2") %>%
#   dplyr::filter(description == "Cancer code, self-reported")
# df4 <- df3 %>%
#   create_pval_category(
#     gene_val = "TSC2",
#     description_val =  "Cancer code, self-reported"
#   )
# df4
# length(df4$pval_category)
# df4 %>%
#   perform_polyphen_score_signif_assoc_by_p_val_category(
#     cancer_type = "Bowel cancer in the colon or rectum",
#     gene_val = "TSC2"
#   )

# t2 <- perform_polyphen_score_signif_assoc_by_p_val_category(
#   df = df21,
#   cancer_type = cancer_type,
#   gene_val = "BRCA1"
# )
# t2
#
# bind_rows(t1, t2) %>%
#   dplyr::mutate(pval_cat = paste0(group1, group2)) %>%
#   dplyr::group_by(pval_cat) %>%
#   dplyr::arrange(pval_cat) %>%
#   dplyr::mutate(adjusted_p_value = p.adjust(p_value, method = "BH"))

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

# unique(significan_test_df$group1)
# df6 <- significan_test_df %>%
#   dplyr::filter(cancer_type == "Breast cancer") %>%
#   dplyr::mutate(pval_cat = paste0(group1, group2)) %>%
#   dplyr::group_by(pval_cat) %>%
#   dplyr::arrange(pval_cat) %>%
#   dplyr::mutate(adjusted_p_value = p.adjust(p_value, method = "BH"))



cancer_types <- unique(significan_test_df$cancer_type)
cancer_types
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
  dplyr::filter(cancer_type == "Breast cancer") 
# %>%
#   dplyr::filter(adjusted_p_value < 0.5)
unique(df2$pval_cat)
df2 %>%
  filter(gene == "BRCA2") %>%
  filter(pval_cat == "[0, 0.01]]0.01, 0.1]")

filename <- "20240225_Polyphen_Genes_Significance_Tests_Cancer_type.csv"
write_csv_data(df = significan_test_adj_pval_df, filename = filename)


title_txt <- "Breast Cancer, Polyphen Score Significance Association Test"
pt <- get_pval_category_adj_pval_scatter_plot(df = df2,
                                        title_txt = title_txt,
                                        max_overlaps_val = 50
                                        )
ppt_filename <- "20240225_Polyphen_Significance_Association_Test.pptx"
save_to_ppt(filename = ppt_filename, pt)

title_txt <- "All Cancers, Polyphen Score Significance Association Test"
pt <- get_pval_category_adj_pval_scatter_plot(df = significan_test_adj_pval_df,
                                              title_txt = title_txt,
                                              max_overlaps_val = 200
)
pt
