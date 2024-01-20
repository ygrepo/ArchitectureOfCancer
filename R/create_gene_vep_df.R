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
source("R/plot_lib.R")
source("R/util_lib.R")

gene <- "APC"

filename <- paste0(gene,"_files_search.txt")
filename <- paste0("genes","/",gene,"/",filename)
file_list <- get_vep_files(filename)

result_list <- lapply(file_list, function(filename) {
  read_vep_file(filename, gene_target)
})
vep_df <- dplyr::bind_rows(result_list)
print(head(vep_df))

filename <- paste0(gene, "_polyphen.csv")
write_csv_data(df = df, filename = filename)

