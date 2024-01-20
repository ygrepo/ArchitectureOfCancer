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

# Get the gene name from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No gene name provided. Usage: Rscript create_gene_vep_df.R <gene_name>", call. = FALSE)
}
gene <- args[1]

filename <- paste0(gene, "_files_search.txt")
filename <- paste0(gene, "/", filename)
file_list <- get_vep_files(filename)

result_list <- lapply(file_list, function(filename) {
  read_vep_file(filename, gene_target)
})
vep_df <- dplyr::bind_rows(result_list)
print(head(vep_df))

filename <- "polyphen.csv"
write_csv_gene_df(df = df, 
                  gene = gene,
                  filename = filename)
