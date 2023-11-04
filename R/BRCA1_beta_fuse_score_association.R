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

cancer_filenames <- getCancerFiles("data/t0.01", gene_target)
data <- as_tibble(read.table(cancer_filenames[["BreastCancer"]],
  header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+"))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))

data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange") %>%
  mutate(ClinVarLabelP = if_else(ClinVarLabel == "LP/P", ProteinChange, "")) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# BRCA1, Breast Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000003-a-2"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Breast Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 5)
ybreaks <- seq(-1.1, 5, by = 1)

pt <- getScatterPlot(study.data,
  title_txt,
  xlabel,
  ylabel,
  x_col,
  y_col,
  alpha = 3,
  point_size = 1,
  point_color = "#00AFBB",
  title_font_size = 8,
  x_y_font_size = 12,
  annotate_text_size = 2,
  annotate.point = TRUE,
  xlimits = xlimits,
  xbreaks = xbreaks,
  ylimits = ylimits,
  ybreaks = ybreaks
)
pt
filename <- "20231104_00000003-a-2_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)


# BRCA1, Breast Cancer, 00000003-b-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000003-b-2"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Breast Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 5)
ybreaks <- seq(-1.1, 5, by = 1)


pt <- getScatterPlot(study.data,
  title_txt,
  xlabel,
  ylabel,
  x_col,
  y_col,
  alpha = 3,
  point_size = 1,
  point_color = "#00AFBB",
  title_font_size = 8,
  x_y_font_size = 12,
  annotate_text_size = 2,
  annotate.point = TRUE,
  xlimits = xlimits,
  xbreaks = xbreaks,
  ylimits = ylimits,
  ybreaks = ybreaks
)
pt
filename <- "20231104_00000003-b-2_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

# BRCA1, Breast Cancer, 00000097-0-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000097-0-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Breast Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 5)
ybreaks <- seq(-1.1, 5, by = 1)


pt <- getScatterPlot(study.data,
  title_txt,
  xlabel,
  ylabel,
  x_col,
  y_col,
  alpha = 3,
  point_size = 1,
  point_color = "#00AFBB",
  title_font_size = 8,
  x_y_font_size = 12,
  annotate_text_size = 2,
  annotate.point = TRUE,
  xlimits = xlimits,
  xbreaks = xbreaks,
  ylimits = ylimits,
  ybreaks = ybreaks
)
pt
filename <- "20231104_00000097-0-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

# BRCA1, Breast Cancer, 00000081-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000081-a-2"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Breast Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 5)
ybreaks <- seq(-1.1, 5, by = 1)


pt <- getScatterPlot(study.data,
  title_txt,
  xlabel,
  ylabel,
  x_col,
  y_col,
  alpha = 3,
  point_size = 1,
  point_color = "#00AFBB",
  title_font_size = 8,
  x_y_font_size = 12,
  annotate_text_size = 2,
  annotate.point = TRUE,
  # xlimits = NULL,
  # xbreaks = NULL,
  # ylimits = NULL,
  # ybreaks = NULL
  xlimits = xlimits,
  xbreaks = xbreaks,
  ylimits = ylimits,
  ybreaks = ybreaks
)
pt
filename <- "20231104_00000081-a-2_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)


# BRCA1, Colorectal Cancer ----

data <- as_tibble(read.table(cancer_filenames[["ColorectalCancer"]],
                            header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+"))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))

data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange") %>%
  mutate(ClinVarLabelP = if_else(ClinVarLabel == "LP/P", ProteinChange, "")) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# BRCA1, Colorectal Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000003-a-2"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Colorectal Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 6)
ybreaks <- seq(-1.1, 6, by = 1)

pt <- getScatterPlot(study.data,
                     title_txt,
                     xlabel,
                     ylabel,
                     x_col,
                     y_col,
                     alpha = 3,
                     point_size = 1,
                     point_color = "#00AFBB",
                     title_font_size = 8,
                     x_y_font_size = 12,
                     annotate_text_size = 2,
                     annotate.point = TRUE,
                     # xlimits = NULL,
                     # xbreaks = NULL,
                     # ylimits = NULL,
                     # ybreaks = NULL
                     
                     xlimits = xlimits,
                     xbreaks = xbreaks,
                     ylimits = ylimits,
                     ybreaks = ybreaks
)
pt
filename <- "20231104_colorectal_00000003-a-2_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)


# BRCA1, Colorectal Cancer, 00000003-b-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000003-b-2"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Colorectal Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 6)
ybreaks <- seq(-1.1, 6, by = 1)


pt <- getScatterPlot(study.data,
                     title_txt,
                     xlabel,
                     ylabel,
                     x_col,
                     y_col,
                     alpha = 3,
                     point_size = 1,
                     point_color = "#00AFBB",
                     title_font_size = 8,
                     x_y_font_size = 12,
                     annotate_text_size = 2,
                     annotate.point = TRUE,
                     xlimits = xlimits,
                     xbreaks = xbreaks,
                     ylimits = ylimits,
                     ybreaks = ybreaks
)
pt
filename <- "20231104_colorectal_00000003-b-2_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

# BRCA1, Colorectal Cancer, 00000097-0-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000097-0-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Colorectal Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 6)
ybreaks <- seq(-1.1, 6, by = 1)

pt <- getScatterPlot(study.data,
                     title_txt,
                     xlabel,
                     ylabel,
                     x_col,
                     y_col,
                     alpha = 3,
                     point_size = 1,
                     point_color = "#00AFBB",
                     title_font_size = 8,
                     x_y_font_size = 12,
                     annotate_text_size = 2,
                     annotate.point = TRUE,
                     xlimits = xlimits,
                     xbreaks = xbreaks,
                     ylimits = ylimits,
                     ybreaks = ybreaks
)
pt
filename <- "20231104_colorectal_00000097-0-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

# BRCA1, Colorectal Cancer, 00000081-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000081-a-2"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "BRCA1, Colorectal Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 3)
xbreaks <- seq(-2, 3, by = 1)
ylimits <- c(-1.1, 6)
ybreaks <- seq(-1.1, 6, by = 1)


pt <- getScatterPlot(study.data,
                     title_txt,
                     xlabel,
                     ylabel,
                     x_col,
                     y_col,
                     alpha = 3,
                     point_size = 1,
                     point_color = "#00AFBB",
                     title_font_size = 8,
                     x_y_font_size = 12,
                     annotate_text_size = 2,
                     annotate.point = TRUE,
                     # xlimits = NULL,
                     # xbreaks = NULL,
                     # ylimits = NULL,
                     # ybreaks = NULL
                     xlimits = xlimits,
                     xbreaks = xbreaks,
                     ylimits = ylimits,
                     ybreaks = ybreaks
)
pt
filename <- "20231104_colorectal_00000081-a-2_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

