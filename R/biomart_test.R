
library(biomaRt)
library(tidyverse)
library(dplyr)
library(kableExtra)
library(export)
library(RColorBrewer)

rm(list = ls())

data_path <- "outputs/data"
figure_path <- "outputs/figures/"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")
source("R/util_lib.R")
source("R/FUSE_lib.R")


mave_data <- readMaveData("data/mave_data_brn_v3.csv")


# install.packages("devtools")
# devtools::install_version("dbplyr", version = "2.3.4")

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes <- c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'sift_score_2076', 'polyphen_score_2076')
peptide_ids <- c('ENSP00000493422')
filters <- 'ensembl_peptide_id'
getBM(attributes = attributes, 
                filters = filters, 
                values = values, 
                mart = ensembl)

gene_ids <- c("ENSG00000130234")  # replace with your gene IDs
transcript_ids <- c("ENST00000252519")  # replace with your transcript IDs
attributes = c('ensembl_gene_id', 'ensembl_peptide_id', 'sift_score_2076', 'polyphen_score_2076')
getBM(attributes = attributes, 
                filters = c('ensembl_gene_id', 'ensembl_transcript_id'), 
                values = list(gene_ids, transcript_ids), 
                mart = ensembl)

gene_ids <- c("ENSG00000130234") 
transcript_ids <- c("ENST00000252519","ENST00000427411")
transcript_ids <- c("ENST00000427411")
attributes = c('ensembl_gene_id', 'ensembl_peptide_id', 'ensembl_transcript_id','sift_score_2076', 'polyphen_score_2076')
getBM(attributes = attributes, 
      filters = c('ensembl_gene_id', 'ensembl_transcript_id'), 
      values = list(gene_ids, transcript_ids), 
      mart = ensembl)

gene_ids <- c("ENSG00000176695") 
transcript_ids <- c("ENST00000618231")
peptide_ids <- c('ENSP00000493422')
attributes = c('ensembl_gene_id', 'ensembl_peptide_id', 'ensembl_transcript_id','sift_score_2076', 'polyphen_score_2076')
getBM(attributes = attributes, 
      filters = c('ensembl_gene_id', 'ensembl_peptide_id','ensembl_transcript_id'), 
      values = list(gene_ids, peptide_ids, transcript_ids), 
      mart = ensembl)

gene_ids <- c("ENSG00000176695") 
strands <- c("ENST00000618231")
attributes = c('ensembl_gene_id', 'strand','sift_score_2076', 'polyphen_score_2076')
getBM(attributes = attributes, 
      filters = c('ensembl_gene_id', 'strand'), 
      values = list(gene_ids, strands), 
      mart = ensembl)


biomaRt::searchAttributes(ensembl)
  attributePages(ensembl)

filters <- listFilters(ensembl)  
attributes <- listAttributes(ensembl)
