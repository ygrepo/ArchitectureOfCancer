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

# MSH2, Colon Cancer ----

gene_target <- "MSH2"

cancer_filenames <- getCancerFiles("data/t0.01", gene_target)
data <- as_tibble(read.table(cancer_filenames[["ColonCancer"]],
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
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
    ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2, Colon Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, Colon Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 9.9)
ybreaks <- seq(-1.1, 9.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_Colon_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")


# MSH2, Colorectal Cancer ----

gene_target <- "MSH2"

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
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
    ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2, Colorectal Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, Colorectal Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 9.9)
ybreaks <- seq(-1.1, 9.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_Colorectal_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")


# MSH2, FHBowel Cancer ----

gene_target <- "MSH2"

data <- as_tibble(read.table(cancer_filenames[["FHBowelCancer"]],
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
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
    ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2, FHBowel Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, FH Bowel Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 7.9)
ybreaks <- seq(-1.1, 7.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_FHBowel_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")


# MSH2, Large Bowel Cancer ----

gene_target <- "MSH2"

data <- as_tibble(read.table(cancer_filenames[["LargeBowelCancer"]],
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
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
    ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2, LargeBowel Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, Large Bowel Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 9.9)
ybreaks <- seq(-1.1, 9.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_LargeBowel_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")


# MSH2, Rectal Cancer ----

gene_target <- "MSH2"

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
  inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange") %>%
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
    ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2, Rectal Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, Rectal Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 7.9)
ybreaks <- seq(-1.1, 7.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_rectal_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")

# MSH2, Small Intestine Cancer ----

gene_target <- "MSH2"

data <- as_tibble(read.table(cancer_filenames[["SmallIntestineCancer"]],
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
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
                                 ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2, Small Intestine Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, Small Intestine Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 9.9)
ybreaks <- seq(-1.1, 9.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_small_intestine_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")


# MSH2,Prostate Cancer ----

gene_target <- "MSH2"

data <- as_tibble(read.table(cancer_filenames[["ProstateCancer"]],
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
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
                                 ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2,Prostate Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, Prostate Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 5.9)
ybreaks <- seq(-1.1, 5.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_prostate_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")


# MSH2, Squamous Cell Cancer ----

gene_target <- "MSH2"

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
  mutate(ClinVarLabelP = if_else(ClinVarLabel %in% c("LB/B", "LP/P"),
                                 ProteinChange, ""
  )) %>%
  mutate(ClinVarLabelP = gsub("p.", "", ClinVarLabelP))


unique(data$id)

# Compute correlation and p-value
studies.data <- getDataByStudy(data)

# MSH2,  Squamous Cell Cancer, 00000003-a-2 ----
# Compute correlation and p-value
study <- "urn:mavedb:00000050-a-1"
study.data <- studies.data$Data[studies.data$Study == study][[1]] %>%
  filter_benign_pathogenic()

res_corr <- get_beta_FUSE_correlation_p_value(study.data)


study_abb <- gsub("urn:mavedb:", "", study)
study_abb
title_txt <- bquote(paste(
  "MSH2, Squamous Cell Cancer, study=", .(study_abb), ", ", beta, " vs. FUSE Score, ", rho, "=",
  .(round(res_corr$corr, 2)), " , p=", .(round(res_corr$pval, 4))
))

xlabel <- "FUSE Score"
ylabel <- "Beta"
x_col <- "FUSE_score"
y_col <- "beta"

xlimits <- c(-3, 2)
xbreaks <- seq(-3, 2, by = 1)
ylimits <- c(-1.1, 0.9)
ybreaks <- seq(-1.1, 0.9, by = 1)

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
  annotate_text_size = 4,
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
pt <- pt + theme(legend.position = "top")
pt
filename <- "20231108_MSH2_squamous_cell_00000050-a-1_beta_fuse.png"
file_path <- paste0(figure_path, filename)
ggsave(file_path, dpi = 600, bg = "white")
