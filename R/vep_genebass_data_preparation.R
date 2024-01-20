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

source("R/io_utils.R")
source("R/plot_lib.R")
source("R/util_lib.R")

gene_target <- "APC"
tissue <- "BreastCancer"
cancer_filenames <- getCancerFiles(paste0("/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/Annotated_Genebass/data/", gene_target), gene_target)
genebass_output <- read_genebass_data_with_variant_process(cancer_filenames = cancer_filenames,
                                                           tissue = tissue)
print(head(genebass_output))


filename <- paste0(gene_target,"_files_search.txt")
filename <- paste0(gene_target,"/",filename)
file_list <- get_vep_files(filename)

result_list <- lapply(file_list, function(filename) {
  read_vep_file(filename, gene_target)
})
vep_df <- dplyr::bind_rows(result_list)
print(head(vep_df))

df <- genebass_output %>%
  inner_join(vep_df, by="variant_id")

filename <- "BRCA1_breast_cancer_polyphen.csv"
write_csv_data(df = df, filename = filename)


