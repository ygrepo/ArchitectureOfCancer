library(dplyr)
library(readr)

figure_path <- "outputs/figures/"
result_data_path <- "outputs/data/"
slides_path <- "outputs/slides/"

readMaveData <- function(file.path_val, gene_target = NULL) {
  df <- as_tibble(read.table(file.path_val,
    header = TRUE, sep = ","
  )) %>%
    mutate(ProteinChange = str_extract(ProteinChange, "(?<=:).+")) %>%
    mutate(
      ClinVarLabel = case_when(
        str_detect(ClinVar.Variant.Category, "US") ~ "US",
        str_detect(ClinVar.Variant.Category, "LB/B") ~ "LB/B",
        str_detect(ClinVar.Variant.Category, "LP/P") ~ "LP/P",
        ClinVar.Variant.Category == "" ~ "None",
        ClinVar.Variant.Category == "-" ~ "None",
        TRUE ~ ClinVar.Variant.Category
      )
    ) %>%
    mutate(ClinVarLabel = trimws(ClinVarLabel)) %>%
    mutate(ClinVarLabel = factor(ClinVarLabel))


  if (!is.null(gene_target)) {
    df <- df %>% filter(gene == gene_target)
  }

  return(df)
}

getCorrelation <- function(label,
                           filename,
                           mave_data,
                           corr_method = "spearman") {
  data <- as_tibble(read.csv(filename))


  data <- data %>%
    mutate(Mutation = str_extract(hgvsp, "(?<=:).+")) %>%
    select(pval, beta, Mutation)

  common_mutations <- intersect(
    unique(data$Mutation),
    unique(mave_data$Mutation)
  )
  cat(paste0(
    "Common mutations:",
    length(unique(common_mutations)), "\n"
  ))

  data_mave_df <- data %>%
    inner_join(mave_data, by = "Mutation") %>%
    select(gene, Mutation, pval, beta, FUSE_score) %>%
    distinct()
  res <- cor.test(data_mave_df$beta, data_mave_df$FUSE_score,
    method = corr_method
  )
  cor_val <- round(res$estimate, 2)[[1]]

  # Create a list with the label and correlation
  result <- list(cancer = label, correlation = cor_val)

  return(result)
}


getCancerFiles <- function(directory_path, cancer.pattern) {
  # List all files in the directory that match the pattern cancer.pattern
  matching_files <- list.files(
    path = directory_path,
    pattern = paste0(cancer.pattern, "_.*Cancer.csv")
  )

  cancer_filenames <- list()

  # Iterate over the matching files and build the mapping
  for (file in matching_files) {
    # Extract the cancer label from the filename
    pattern.toextract <- paste0(cancer.pattern, "_(.*).csv")
    label <- sub(pattern.toextract, "\\1", file)
    # label <- sub("BRCA1_(.*).csv", "\\1", file)


    # Create the full file path
    file_path <- file.path(directory_path, file)

    # Add the label and file path to the mapping
    cancer_filenames[[label]] <- file_path
  }

  cancer_filenames
}

get_vep_files <- function(directory_path, file_pattern) {
  # List all files in the directory that match the pattern cancer.pattern
  matching_files <- list.files(
    path = directory_path,
    pattern = file_pattern
  )

  cancer_filenames <- list()

  # Iterate over the matching files and build the mapping
  for (file in matching_files) {
    # Create the full file path
    file_path <- file.path(directory_path, file)
    # Add the label and file path to the mapping
    cancer_filenames <- c(cancer_filenames, file_path)
  }

  return(cancer_filenames)
}


read_vep_file <- function(filename, gene) {
  vep_output <- as_tibble(read.table(filename,
    header = TRUE, sep = "\t"
  )) %>%
    dplyr::filter((gene == gene)) %>%
    # dplyr::filter((gene == gene_target) &
    #                 (consequence %in% unique(genebass_output$genebass_consequence))) %>%
    # filter(substr(markerID, nchar(markerID) - 2, nchar(markerID)) == "T/C") %>%
    dplyr::mutate(markerID = substr(markerID, 7, nchar(markerID))) %>%
    dplyr::mutate(markerID = stringr::str_replace_all(markerID, "_([A-Za-z])\\/([A-Za-z])", "-\\1-\\2")) %>%
    # mutate(markerID = str_replace(markerID, "_T/C$", "-T-C"))%>%
    dplyr::rename("variant_id" = "markerID") %>%
    dplyr::rename("vep_consequence" = "consequence") %>%
    dplyr::arrange(variant_id) %>%
    dplyr::mutate(
      polyphen_label = stringr::str_extract(polyphen, "^[a-zA-Z]+"), # Extracting letters before the parentheses
      polyphen_score = as.numeric(stringr::str_extract(polyphen, "\\d+\\.\\d+")) # Extracting floating point number
    ) %>%
    dplyr::mutate(
      polyphen_label_simplified = case_when(
        str_detect(polyphen_label, "possibly") ~ "probably",
        is.na(polyphen_label) ~ "unknown",
        TRUE ~ polyphen_label
      )
    ) %>%
    dplyr::rename("polyphen_pval" = "Pvalue")

  return(vep_output)
}

read_genebass_data_with_variant_process <- function(tissue, cancer_filenames) {
  genebass_output <- as_tibble(read.table(cancer_filenames[[tissue]],
    header = TRUE, sep = ","
  )) %>%
    # filter(consequence == "missense_variant") %>%
    # filter(substr(variant_id, nchar(variant_id) - 2, nchar(variant_id)) == "T-C") %>%
    dplyr::mutate(variant_id = substr(variant_id, 4, nchar(variant_id))) %>%
    dplyr::rename("genebass_consequence" = "consequence") %>%
    dplyr::arrange(variant_id) %>%
    dplyr::rename("genebass_pval" = "pval")

  return(genebass_output)
}


write_csv_data <- function(df, filename, delim = " ") {
  setwd("~/github/ArchitectureOfCancer/")
  filename <- paste0(result_data_path, filename)
  print(paste0("Saving data to:", filename))
  readr::write_delim(df, filename, delim = delim)
}

read_csv_data <- function(filename, delim = " ") {
  setwd("~/github/ArchitectureOfCancer")
  filename <- paste0(result_data_path, filename)
  return(dplyr::as_tibble(readr::read_delim(filename, delim = delim)))
}


print_html_df <- function(df,
                          caption_txt,
                          file_path = NULL) {
  rep <- df %>%
    kbl(caption = caption_txt) %>%
    kable_classic(
      bootstrap_options = c("stripped", "condensed", "responsive"),
      full_width = F,
      html_font = "Arial"
    )
  if (!is.null(file_path)) {
    rep %>% save_kable(file = file_path, self_contained = TRUE)
  }
  return(rep)
}


save_to_ppt <- function(filename, pt) {
  setwd("~/github/ArchitectureOfCancer")
  filename <- paste0(slides_path, filename)
  print(paste0("Saving to:", filename))
  graph2ppt(pt, filename, width = 7, height = 7, append = TRUE)
  
}
