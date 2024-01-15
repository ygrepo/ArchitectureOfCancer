library(tidyverse)
library(dplyr)
library(ggcorrplot)
library(plotly)
library(kableExtra)
library(export)
library(RColorBrewer)
library(ggpubr)

rm(list = ls())

data_path <- "outputs/data"
figure_path <- "outputs/figures/"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")
source("R/util_lib.R")

gene_target <- "BRCA1"
tissue <- "BreastCancer"
cancer_filenames <- getCancerFiles(paste0("data/", gene_target), gene_target)
genebass_output <- read_genebass_data_with_variant_process(cancer_filenames = cancer_filenames,
                                                           tissue = tissue)
file_list <- get_vep_files(directory_path = "data/vep_outputs/breast_cancer", 
                           file_pattern = "cancers")

result_list <- lapply(file_list, function(filename) {
  read_vep_file(filename, gene_target)
})
vep_df <- dplyr::bind_rows(result_list)

df <- genebass_output %>%
  inner_join(vep_df, by="variant_id")

filename <- "BRCA1_breast_cancer_polyphen.csv"
write_csv_data(df = df, filename = filename)


