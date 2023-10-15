rm(list = ls())

library(tidyverse)
library(plotly)
library(kableExtra)
library(data.table)
library(export)

setwd("~/github/ArchitectureOfCancer")

ouput_path <- "ouputs"
figure_path <- "outputs/figures/"

setwd("~/github/ArchitectureOfCancer/")

source("R/plot_lib.R")
source("R/util_lib.R")



dt <- fread("data/gene-phewas-exomes_ENSG00000095002_MSH2.csv")
filter_condition <- function(x) grepl("cancer|carcinoma|melanoma", 
                                      x, ignore.case = TRUE)
filtered_dt <- dt[filter_condition(Description)]

# Load MAVE data ----

gene_target <- "MSH2"
mave_data <- as_tibble(read.table("data/mave_data_brn_v2.csv",
  header = TRUE, sep = ","
)) %>%
  filter(gene == gene_target)

head(mave_data)
glimpse(mave_data)
mave_data <- mave_data %>%
  mutate(Mutation = str_extract(ProteinChange, "(?<=:).+"))
unique(mave_data$id)

# All being equals mutations can have different FUSE score
# df <- mave_data %>%
#   filter(Mutation == "p.Pro34Leu")

mave_data <- mave_data %>% select(id, gene, FUSE_score, Mutation)
studies <- unique(mave_data$id)

cancer_filenames <- getCancerFiles("data/t0.01", "MSH2")

corr_dt <- data.table(
  Study = character(0), # Numeric column (replace with your desired data type)
  Cancer = character(0), # Character column (replace with your desired data type)
  Cor = numeric(0) # Logical column (replace with your desired data type)
)


for (study in studies) {
  print(study)
  df <- mave_data %>% filter(id == study)
  print(head(df))

  # Create an empty data table to store the results for the current study
  study_results <- data.table(
    Study = character(0),
    Cancer = character(0),
    Cor = numeric(0)
  )

  for (label in names(cancer_filenames)) {
    filename <- cancer_filenames[[label]]
    result <- getCorrelation(label, filename, df)
    study_results <- rbind(
      study_results,
      data.table(
        Study = study,
        Cancer = label,
        Cor = result$correlation
      )
    )
  }

  # Append the study_results to corr_dt
  corr_dt <- rbind(corr_dt, study_results)

}


setorder(corr_dt, Cancer)
corr_dt$Study <-
  gsub("^urn:mavedb:", "", corr_dt$Study)

glimpse(corr_dt)

# Set scipen to a large value to avoid scientific notation
# options(scipen = 999)

font_size <- 10
legend_font_size <- 8
xlimits <- c(-0.4, 0.4)
xbreaks <- round(seq(-0.4, 0.4, by = 0.1), digits = 2)
expression_val <- expression("MSH2: Correlations between " ~ beta ~
  " and FUSE Score by Cancer Type (" ~ beta ~ " p-val < 0.01)")
pt <- getCorrelationByCancerTypePlot(corr_dt,
  lab_x_txt = "Spearman Correlation",
  expression_val,
  xlimits = xlimits,
  xbreaks = xbreaks,
  point_size = 2,
  point_alpha = 0.8
)
pt

filename <- "20231015_MSH2_BETA_LESS_0.01_FUSE_Correlations_Cancer.png"
file_path <- paste0(figure_path, filename)
png(file_path, width = 800, height = 600, units = "px", pointsize = 12)
print(pt)
dev.off()
# graph2png(x = pt, file = file_path, dpi = 600, aspectr = 2.2)


