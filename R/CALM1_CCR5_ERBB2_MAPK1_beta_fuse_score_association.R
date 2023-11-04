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

# CALM1, RCC Cancer ----

gene_target <- "CALM1"

cancer_filenames <- getCancerFiles("data/t0.01", gene_target)
data <- as_tibble(read.table(cancer_filenames[["KidneyRenalCellCancer"]],
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

# CALM1, RCC Cancer, 00000001-c-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000001-c-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "CALM1, RCC Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 11)
ybreaks <- seq(-1.1, 11, by = 1)

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
  # xlimits = xlimits,
  # xbreaks = xbreaks,
  # ylimits = ylimits,
  # ybreaks = ybreaks
  xlimits = NULL,
  xbreaks = NULL,
  ylimits = NULL,
  ybreaks = NULL
  
)
pt
filename <- "20231104_CALM1_RCC_00000001-c-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)


# CCR5, FH Breast Cancer ----

gene_target <- "CCR5"
cancer_filenames <- getCancerFiles("data/t0.01", gene_target)
data <- as_tibble(read.table(cancer_filenames[["FHBreastCancer"]],
                             header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+"))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))

data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), 
             by = "ProteinChange", relationship = "many-to-many") %>%
  mutate(ClinVarLabelP = if_else(ClinVarLabel == "LP/P", ProteinChange, "")) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# CCR5, FH Breast Cancer, 00000047-a-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000047-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "CCR5, FH Breast Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 1)
xbreaks <- seq(-2, 1, by = 1)
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
                     # xlimits = NULL,
                     # xbreaks = NULL,
                     # ylimits = NULL,
                     # ybreaks = NULL
                     
)
pt
filename <- "20231104_CCR5_FHBC_00000047-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

# CCR5, FH Breast Cancer, 00000047-b-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000047-b-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "CCR5, FH Breast Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 1)
xbreaks <- seq(-2, 1, by = 1)
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
                     # xlimits = NULL,
                     # xbreaks = NULL,
                     # ylimits = NULL,
                     # ybreaks = NULL
                     
)
pt
filename <- "20231104_CCR5_FHBC_00000047-b-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

# CCR5, FH Breast Cancer, 00000047-c-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000047-c-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "CCR5, FH Breast Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 1)
xbreaks <- seq(-2, 1, by = 1)
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
                     # xlimits = NULL,
                     # xbreaks = NULL,
                     # ylimits = NULL,
                     # ybreaks = NULL
                     
)
pt
filename <- "20231104_CCR5_FHBC_00000047-c-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)

# ERBB2, Rectal Cancer ----

gene_target <- "ERBB2"
cancer_filenames <- getCancerFiles("data/t0.01", gene_target)
data <- as_tibble(read.table(cancer_filenames[["RectalCancer"]],
                             header = TRUE, sep = ","
)) %>%
  mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+"))


common_variants <- intersect(
  unique(data$ProteinChange),
  unique(mave_data$ProteinChange)
)

print(length(unique(common_variants)))

data <- data %>%
  inner_join(mave_data %>% filter(gene == gene_target), 
             by = "ProteinChange", relationship = "many-to-many") %>%
  mutate(ClinVarLabelP = if_else(ClinVarLabel == "LP/P", ProteinChange, "")) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# ERBB2, Rectal Cancer, 00000051-b-1" ----
# Compute correlation and p-value
study <- "urn:mavedb:00000051-b-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "ERBB2, Rectal Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-2, 1)
xbreaks <- seq(-2, 1, by = 1)
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
                     # xlimits = xlimits,
                     # xbreaks = xbreaks,
                     # ylimits = ylimits,
                     # ybreaks = ybreaks
                     xlimits = NULL,
                     xbreaks = NULL,
                     ylimits = NULL,
                     ybreaks = NULL
                     
)
pt
filename <- "20231104_ERBB2_rectal_00000051-b-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)


# MAPK1, C43 Malignant Skin Cancer ----

gene_target <- "MAPK1"

cancer_filenames <- getCancerFiles("data/t0.01", gene_target)
data <- as_tibble(read.table(cancer_filenames[["C43MalignantSkinCancer"]],
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

# MAPK1, C43 Malignant Skin Cancer, 00000103-b-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000103-b-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MAPK1, Skin Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 11)
ybreaks <- seq(-1.1, 11, by = 1)

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
                     # xlimits = xlimits,
                     # xbreaks = xbreaks,
                     # ylimits = ylimits,
                     # ybreaks = ybreaks
                     xlimits = NULL,
                     xbreaks = NULL,
                     ylimits = NULL,
                     ybreaks = NULL
                     
)
pt
filename <- "20231104_MAPK1_skin_00000103-b-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)


# MAPK1, Squamous Cell Cancer ----

gene_target <- "MAPK1"

data <- as_tibble(read.table(cancer_filenames[["SquamousCellCancer"]],
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

# MAPK1, Squamous CellCancer, 00000103-b-1 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000103-b-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]]
res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MAPK1, Squamous Cell Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 11)
ybreaks <- seq(-1.1, 11, by = 1)

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
                     # xlimits = xlimits,
                     # xbreaks = xbreaks,
                     # ylimits = ylimits,
                     # ybreaks = ybreaks
                     xlimits = NULL,
                     xbreaks = NULL,
                     ylimits = NULL,
                     ybreaks = NULL
                     
)
pt
filename <- "20231104_MAPK1_squamous_cell_00000103-b-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
graph2png(pt, file = file_path, dpi = 600, aspectr = 1.2)
