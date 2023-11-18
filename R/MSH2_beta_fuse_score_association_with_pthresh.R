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

gene_target <- "MSH2"
pval_threshold <- 0.05

cancer_filenames <- getCancerFiles("data/MSH2", gene_target)

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
  "MSH2, Breast Cancer, ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 2)
xbreaks <- seq(-2, 2, by = 1)
ylimits <- c(-1, 5)
ybreaks <- seq(-1, 5, by = 1)


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
filename <- "20231118_MSH2_Beta_Fuse.pptx"
file_path <- paste0(figure_path, filename)
graph2ppt(pt, file_path, width = 7, height = 6, append = TRUE)
# ggsave(file_path, dpi = 600, bg="white")


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
  "MSH2, Colorectal Cancer, ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 1)
xbreaks <- seq(-3, 1, by = 1)
ylimits <- c(2, 8)
ybreaks <- seq(3, 8, by = 1)


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

graph2ppt(pt, file_path, width = 7, height = 5, append = TRUE)


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
  "MSH2, Lung Cancer, ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 2)
xbreaks <- seq(-2, 2, by = 1)
ylimits <- c(2, 7)
ybreaks <- seq(2, 7, by = 1)


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

graph2ppt(pt, file_path, width = 7, height = 5, append = TRUE)



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



data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange")

res_corr <- get_beta_FUSE_correlation_p_value(data)
title_txt <- bquote(paste(
  "MSH2, Prostate Cancer, ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-1, 2)
xbreaks <- seq(-1, 2, by = 1)
ylimits <- c(2, 6)
ybreaks <- seq(2, 6, by = 1)


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

graph2ppt(pt, file_path, width = 7, height = 5, append = TRUE)


