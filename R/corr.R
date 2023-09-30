rm(list = ls())
# install.packages("kableExtra")
# install.packages("tidyverse")
# BiocManager::install("biomaRt")
# library(biomaRt)
# library(VariantAnnotation)
library(tidyverse)
library(ggcorrplot)
library(plotly)
library(kableExtra)

setwd("~/github/ArchitectureOfCancer")

ouput_path <- "ouputs"
figure_path <- "figures"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")

data <- as_tibble(read.table("data/genbass_breastcancer_custom_BRCA1.csv",
  header = TRUE, sep = ","
))

glimpse(data)
data <- data %>%
  mutate(Protein = str_extract(hgvsp, "(?<=:).+"))
unique(data$Protein)

mave_data <- as_tibble(read.table("data/mave_data_brn_v2.csv",
  header = TRUE, sep = ","
))
mave_data <- mave_data %>% filter(gene == "BRCA1")
head(mave_data)
glimpse(mave_data)
mave_data <- mave_data %>%
  mutate(Protein = str_extract(ProteinChange, "(?<=:).+"))
unique(mave_data$Protein)

mave_data <- mave_data %>% select(-1)

common_proteins <- intersect(unique(data$Protein), unique(mave_data$Protein))

data_mave_df <- data %>% inner_join(mave_data, by = "Protein")
cor_data <- data_mave_df %>% select(beta, raw_score)
corrdata <- cor(cor_data, use = "na.or.complete")
p.mat <- cor_pmat(cor_data)

ggcorrplot(corrdata,
  title = "Correlation matrix for Beta And Raw Score",
  lab = TRUE, p.mat = p.mat, sig.level = .05
)

title_txt <- "Mave Raw Score vs. Beta"
xlabel <- "Raw Score"
ylabel <- "Beta"
x_col <- "raw_score"
y_col <- "beta"
xlimits <- c(-3, 1)
xbreaks <- seq(-3, 1, by = 1)
ylimits <- c(-1, 5)
ybreaks <- seq(-1, 5, by = 1)
pt <- get_scatter_plot(data_mave_df,
  title_txt,
  xlabel,
  ylabel,
  x_corr = -1.7,
  y_corr = -1,
  alpha = 2,
  raw_score,
  beta,
  point_size = 2,
  point_color = "#00AFBB",
  font_size = 14,
  xlimits = xlimits,
  xbreaks = xbreaks,
  ylimits = ylimits,
  ybreaks = ybreaks
)
pt

cor_data <- data_mave_df %>% select(beta, FUSE_score)
corrdata <- cor(cor_data, use = "na.or.complete")
p.mat <- cor_pmat(cor_data)
title_txt <- "Mave Raw Score vs. Beta"
xlabel <- "Raw Score"
ylabel <- "Beta"
x_col <- "raw_score"
y_col <- "beta"
xlimits <- c(-3, 1)
xbreaks <- seq(-3, 1, by = 1)
ylimits <- c(-1, 5)
ybreaks <- seq(-1, 5, by = 1)
title_txt <- "Mave Raw Score vs. Beta"
xlabel <- "Raw Score"
ylabel <- "Beta"
x_col <- "raw_score"
y_col <- "beta"
xlimits <- c(-3, 1)
xbreaks <- seq(-3, 1, by = 1)
ylimits <- c(-1, 5)
ybreaks <- seq(-1, 5, by = 1)
pt <- get_scatter_plot(data_mave_df,
                       title_txt,
                       xlabel,
                       ylabel,
                       x_corr = -1.7,
                       y_corr = -1,
                       alpha = 2,
                       raw_score,
                       beta,
                       point_size = 2,
                       point_color = "#00AFBB",
                       font_size = 14,
                       annotate_text_size = 4,
                       xlimits = xlimits,
                       xbreaks = xbreaks,
                       ylimits = ylimits,
                       ybreaks = ybreaks
)
pt
