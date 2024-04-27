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

tissue <- "BreastCancer"
ppt_filename <- "20240117_BRCA1_Polyphen_PVal.pptx"
size_col <- "beta"
xlabel <- "P Value Category"
ylabel <- "Polyphen Score"


# Breast Cancer ----
cancer_type <- "Breast cancer"
gene_val <- "BRCA1"
signifPVal <- (-log10(0.05))
df_cancer <- read_csv_polyphen_df(gene) %>%
  dplyr::filter((description == cancer_type) & (gene == gene_val)) %>%
  na.omit((polyphen_score)) %>%
  na.omit(Pvalue) %>%
  dplyr::mutate(Log10PVal = -log10(Pvalue)) %>%
  dplyr::mutate(variantIdSign = case_when(
    Log10PVal <= signifPVal ~ "",
    TRUE ~ variant_id
  ))

title_txt <- "Breast Cancer, BRCA1, Polyphen Score vs. -LOG10PValue"
font_size <- 12
legend_bottom <- NULL
alpha <- 0.6
pt <- df_cancer %>%
  ggplot2::ggplot(aes(
    x = polyphen_score,
    y = Log10PVal,
  )) +
  geom_point(
    data = df_cancer,
    alpha = alpha,
    shape = 21,
    size = 5,
    fill = "steelblue",
    colour = "black"
  ) +
  labs(
    x = "Polyphen Score",
    y = "-log10(PValue)",
    title = title_txt
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_publication(
    base_size = font_size,
    legend_bottom = legend_bottom
  ) +
  theme(
    plot.title = element_text(
      color = "black",
      size = font_size,
      face = "bold", hjust = 0.5
    ),
    axis.title.x = element_text(size = font_size, face = "bold"),
    axis.title.y = element_text(size = font_size, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_x_continuous(
    breaks = c(seq(0, 1.1, 0.2)),
    limits = c(0, 1.1)
  ) +
  scale_y_continuous(
    breaks = c(seq(0, 2.1, 0.2)),
    limits = c(0, 2.1)
  )

var_label_col <- "variantIdSign"
pt <- pt + geom_text_repel(
  data = df_cancer,
  aes(
    x = polyphen_score,
    y = Log10PVal,
    label = .data[[var_label_col]]
  ),
  size = 3,
  max.overlaps = 50
)
pt
ppt_filename <- "20240427_BreastCancer_BRCA1_Polyphen_LOG10PVal.pptx"

save_to_ppt(filename = ppt_filename, pt)
