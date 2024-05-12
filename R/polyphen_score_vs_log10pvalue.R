library(tidyverse)
library(dplyr)
library(ggcorrplot)
library(ggrepel)
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
tissue <- "BreastCancer"

cancer_filenames <- getCancerFiles(gene_val)
genebass_df <- read_genebass_data_with_variant_process(cancer_filenames = cancer_filenames,
                                                           tissue = tissue)

signifPVal <- (-log10(0.05))
df_cancer <- read_csv_polyphen_df(gene_val) %>%
  dplyr::filter((description == cancer_type) & (gene == gene_val)) %>%
  na.omit((polyphen_score)) %>%
  na.omit(Pvalue) %>%
  dplyr::mutate(Log10PVal = -log10(Pvalue))

df_cancer$AF_bin <- factor(df_cancer$AF_bin, levels = c("other", "rare", "common"))

df_cancer <- df_cancer %>%
  dplyr::left_join(genebass_df) %>%
  dplyr::mutate(variantIdSign = case_when(
    Log10PVal <= signifPVal ~ "",
    TRUE ~ shgvsp
  ))

title_txt <- "Breast Cancer, BRCA1, Polyphen Score vs. -LOG10PValue"
font_size <- 12
legend_bottom <- NULL
alpha <- 1
pt <- df_cancer %>%
  ggplot2::ggplot(aes(
    x = polyphen_score,
    y = Log10PVal,
  )) +
  geom_point(
    aes(
      color = BETA,
      size = AF_bin
    ),
    alpha = alpha,
    shape = 16
  ) +
  labs(
    x = "Polyphen Score",
    y = "-log10(PValue)",
    size = "AF Bin",
    title = title_txt
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  cowplot::theme_cowplot() +
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
  ) +
  scale_size_manual(
    values = c("other" = 3, "rare" = 5, "common" = 10), # Adjust sizes as needed
    name = "AF"
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint = 0,
    name = "BETA"
  )

var_label_col <- "variantIdSign"
pt <- pt + geom_text_repel(
  data = df_cancer,
  aes(
    x = polyphen_score,
    y = Log10PVal,
    label = .data[[var_label_col]]
  ),
  box.padding = .6,
  label.size = 6,
  max.overlaps = 10
)
pt

filename <- "20240509_BreastCancer_BRCA1_Polyphen_LOG10PVal.jpg"
save_plot(filename = filename,
          pt = pt)


save_to_ppt(filename = ppt_filename, pt)
