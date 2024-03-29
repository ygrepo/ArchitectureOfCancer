rm(list = ls())
library(tidyverse)
library(ggcorrplot)
library(plotly)
library(kableExtra)
library(export)


data_path <- "outputs/data"
figure_path <- "outputs/figures/"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")
source("R/util_lib.R")

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

mave_data %>%
  group_by(ClinVar.ClinicalSignificance, ClinVar.Variant.Category, ClinVarLabel) %>%
  select(ClinVar.ClinicalSignificance, ClinVar.Variant.Category, ClinVarLabel) %>%
  unique() %>%
  print_html_df(caption_txt = "Clinvar Significance, Categories And Label")

gene_target <- "BRCA1"
pval_threshold <- 0.05

cancer_filenames <- getCancerFiles("data/BRCA1", gene_target)

# Breast Cancer ----
data <- as_tibble(read.table(cancer_filenames[["BreastCancer"]],
  header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+")) %>%
  filter(pval < pval_threshold) %>%
  mutate(LOG10AF = -log10(allele_frequency)) %>%
  # mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
  #                                ProteinChange, "")) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ProteinChange))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))

data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange")

res_corr <- get_beta_FUSE_correlation_p_value(data)
title_txt <- bquote(paste(
  "BRCA1, Breast Cancer, ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 1)
xbreaks <- seq(-1.5, 1, by = 1)
ylimits <- c(3, 5)
ybreaks <- seq(3, 5, by = 1)


pt <- getScatterPlot(data,
  title_txt,
  xlabel,
  ylabel,
  x_col,
  y_col,
  alpha = 3,
  point_size = 2,
  point_color = "#00AFBB",
  title_font_size = 10,
  x_y_font_size = 12,
  annotate_text_size = 4,
  annotate.point = TRUE,
  xlimits = xlimits,
  xbreaks = xbreaks,
  ylimits = ylimits,
  ybreaks = ybreaks
)
pt
filename <- "20231118_BRCA1_BreastCancer_Beta_Fuse.pptx"
# filename <- "20231118_BRCA1_BreastCancer_Beta_Fuse.png"
file_path <- paste0(figure_path, filename)
graph2office(pt, file_path, width = 6, height = 5)
# ggsave(file_path, dpi = 600, bg="white")
# graph2png(pt, file = file_path, dpi = 600, aspectr = 1)


# Colorectal Cancer ----
data <- as_tibble(read.table(cancer_filenames[["ColorectalCancer"]],
  header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+")) %>%
  filter(pval < pval_threshold) %>%
  mutate(LOG10AF = -log10(allele_frequency)) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ProteinChange))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))

data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange")

res_corr <- get_beta_FUSE_correlation_p_value(data)
title_txt <- bquote(paste(
  "BRCA1, Colorectal Cancer, ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-0.5, 1)
xbreaks <- seq(-0.5, 1, by = 1)
ylimits <- c(5.5, 7.5)
ybreaks <- seq(5.5, 7.5, by = 1)


pt <- getScatterPlot(data,
                     title_txt,
                     xlabel,
                     ylabel,
                     x_col,
                     y_col,
                     alpha = 3,
                     point_size = 2,
                     point_color = "#00AFBB",
                     title_font_size = 10,
                     x_y_font_size = 12,
                     annotate_text_size = 4,
                     annotate.point = TRUE,
                     xlimits = xlimits,
                     xbreaks = xbreaks,
                     ylimits = ylimits,
                     ybreaks = ybreaks
)
pt

filename <- "20231118_BRCA1_ColorectalCancer_Beta_Fuse.pptx"
file_path <- paste0(figure_path, filename)
graph2office(pt, file_path, width = 6, height = 5)


# Lung Cancer ----
data <- as_tibble(read.table(cancer_filenames[["LungCancer"]],
                             header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+")) %>%
  filter(pval < pval_threshold) %>%
  mutate(LOG10AF = -log10(allele_frequency)) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ProteinChange))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))

data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange")

res_corr <- get_beta_FUSE_correlation_p_value(data)
title_txt <- bquote(paste(
  "BRCA1, Lung Cancer, ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 2)
xbreaks <- seq(-0.5, 1, by = 1)
ylimits <- c(2, 8)
ybreaks <- seq(5.5, 7.5, by = 1)


pt <- getScatterPlot(data,
                     title_txt,
                     xlabel,
                     ylabel,
                     x_col,
                     y_col,
                     alpha = 3,
                     point_size = 2,
                     point_color = "#00AFBB",
                     title_font_size = 10,
                     x_y_font_size = 12,
                     annotate_text_size = 4,
                     annotate.point = TRUE,
                     xlimits = xlimits,
                     # xbreaks = xbreaks,
                     ylimits = ylimits
                     # ybreaks = ybreaks
)
pt

filename <- "20231118_BRCA1_LungCancer_Beta_Fuse.pptx"


# Prostate Cancer ----
data <- as_tibble(read.table(cancer_filenames[["ProstateCancer"]],
                             header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+")) %>%
  filter(pval < pval_threshold) %>%
  mutate(LOG10AF = -log10(allele_frequency)) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ProteinChange))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))



