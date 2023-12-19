library(tidyverse)
library(dplyr)
library(ggcorrplot)
library(plotly)
library(kableExtra)
library(export)
library(RColorBrewer)
library(ggpubr)
library(vroom)
rm(list = ls())

data_path <- "outputs/data"
figure_path <- "outputs/figures/"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")
source("R/util_lib.R")
source("R/FUSE_lib.R")


mave_data <- readMaveData("data/mave_data_brn_v3.csv")

df <- mave_data %>%
  dplyr::filter((gene == "BRCA1") & (!is.na(raw_score))) %>%
  dplyr::select(
    gene, aapos, aaref, aaalt,
    raw_score, norm_raw_score,
    pos_score, FUSE_score
  ) %>%
  dplyr::distinct()

filename <- "20231218_mave_BRCA1_for_scoring.csv"
write_csv_data(
  df = df %>%
    dplyr::select(
      gene, aapos, aaref, aaalt,
      raw_score, norm_raw_score
    ) %>%
    dplyr::distinct(),
  filename = filename,
  delim = ","
)

res_df <- get_FUSE_score_df(
  filename = filename,
  pos_mean_method = "js"
)
res_df %>%
  dplyr::filter((aapos == 1) & (aaref == "M") & (aaalt %in% c(
    "R",
    "K",
    "T",
    "V",
    "I"
  ))) %>%
  dplyr::arrange(aaalt)


df %>%
  dplyr::filter((aapos == 1) & (aaref == "M") & (aaalt %in% c(
    "R",
    "K",
    "T",
    "V",
    "I"
  ))) %>%
  dplyr::arrange(aaalt)
