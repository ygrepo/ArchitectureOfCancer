library(tidyverse)
library(ggcorrplot)
library(plotly)
library(kableExtra)
library(export)

rm(list = ls())

data_path <- "outputs/data"
figure_path <- "outputs/figures/"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")
source("R/util_lib.R")
source("R/io_utils.R")

mave_data <- readMaveData("data/mave_data_brn_v3.csv") %>%
  mutate(
    ClinVarLabel = case_when(
      str_detect(ClinVar.Variant.Category, "US") ~ "US",
      str_detect(ClinVar.Variant.Category, "LB/B") ~ "LB/B",
      str_detect(ClinVar.Variant.Category, "LP/P") ~ "LP/P",
      TRUE ~ ClinVar.Variant.Category
    )
  )

mave_data$ClinVarLabel <- trimws(mave_data$ClinVarLabel)
mave_data$ClinVarLabel <- factor(mave_data$ClinVarLabel)


gene_target <- "BRCA1"
pval_threshold <- 0.01

cancer_filenames <- getCancerFiles(gene_target)


# Breast Cancer ----
data <- dplyr::as_tibble(read.table(cancer_filenames[["BreastCancer"]],
  header = TRUE, sep = ","
)) %>%
  dplyr::mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+")) %>%
  # filter(pval < pval_threshold) %>%
  dplyr::mutate(adj_p_val = p.adjust(pval, method = "fdr")) %>%
  dplyr::mutate(LOG10AF = -log10(allele_frequency)) %>%
  dplyr::mutate(LOG10PVAL = -log10(pval)) %>%
  dplyr::mutate(LOG10ADJPVAL = -log10(adj_p_val)) %>%
  #dplyr::filter(LOG10ADJPVAL > 1) %>%
  dplyr::mutate(ClinVarLabelP = gsub("p.", "", ProteinChange)) %>%
  dplyr::mutate(beta_type = case_when(
    beta >= 2 & LOG10ADJPVAL >= -log10(0.05) ~ "up",
    beta <= 1 & LOG10ADJPVAL >= -log10(0.05) ~ "down",
    TRUE ~ "ns"
  )) %>%
  na.omit(LOG10ADJPVAL) %>%
  na.omit(beta)

pt <- get_beta_pval_violin_plot(df = data,
                                title_tx = "BRCA1 Cancer Variants -log10(adj_pval) vs. Beta")
pt
filename <- "20240204_BRCA1_Beta_AdjPVal.pptx"
file_path <- paste0(figure_path, filename)
save_to_ppt(filename = filename,
            pt = pt)

data %>%
  na.omit(beta) %>%
  pull(beta) %>%
  min() %>%
  floor()

df <- data %>%
  filter(beta_type == "up")
# -2

data %>%
  na.omit(beta) %>%
  pull(beta) %>%
  max() %>%
  ceiling()

