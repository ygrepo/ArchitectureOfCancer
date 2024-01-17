library(dplyr)
library(readr)

result_data_path <- "outputs/data/"

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

get_beta_FUSE_correlation_p_value <- function(df, method_txt = "spearman") {
  cor_test <- cor.test(df$beta, df$FUSE_score, method = method_txt)
  cor_val <- cor_test$estimate
  p_val <- cor_test$p.value
  return(list(corr = cor_val, pval = p_val))
}


getDataByStudy <- function(df) {
  result <- tibble(Study = character(), Data = list())

  for (study in unique(df$id)) {
    # Extract ID and Data for the current iteration
    current_data <- df %>% filter(id == study)

    # Create a new row with ID and Data
    new_row <- tibble(Study = as.character(study), Data = list(current_data))

    # Bind the new row to the result tibble
    result <- bind_rows(result, new_row)
  }

  return(result)
}


filter_benign_pathogenic <- function(df) {
  df <- df %>%
    filter(ClinVarLabel %in% c("LB/B", "LP/P"))
  return(df)
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



get_genebass_mave_data_with_pval <- function(gene_target,
                                             mave_data,
                                             cancer_label,
                                             cancer_filenames,
                                             n_outliers = 5) {
  data <- as_tibble(read.table(cancer_filenames[[cancer_label]],
    header = TRUE, sep = ","
  )) %>%
    mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+")) %>%
    mutate(pval_category = case_when(
      pval <= 0.01 ~ "[0, 0.01]",
      (0.01 < pval) & (pval < 0.1) ~ "[0.01, 0.1[",
      (0.1 < pval) ~ "[0.1, 1[",
      TRUE ~ NA_character_
    )) %>%
    mutate(pval_category = factor(pval_category)) %>%
    mutate(LOG10AF = -log10(allele_frequency)) %>%
    mutate(LOG10PVAL = -log10(pval)) %>%
    mutate(ClinVarLabelP = gsub("p.", "", ProteinChange))


  common_variants <- intersect(
    unique(data$ProteinChange),
    unique(mave_data$ProteinChange)
  )

  print(paste0("Common Variants:", length(unique(common_variants))))

  data <- data %>%
    inner_join(mave_data %>% filter(gene == gene_target), by = "ProteinChange")

  data <- data %>%
    filter(!is.na(pval_category)) %>%
    group_by(pval_category) %>%
    mutate(
      lower_bound = quantile(FUSE_score, 0.25) - 1.5 * IQR(FUSE_score),
      upper_bound = quantile(FUSE_score, 0.75) + 1.5 * IQR(FUSE_score),
      is_lower_outlier = ifelse(FUSE_score < lower_bound, TRUE, FALSE),
      is_upper_outlier = ifelse(FUSE_score > upper_bound, TRUE, FALSE),
      is_outlier = ifelse(is_lower_outlier | is_upper_outlier,
        TRUE, FALSE
      ),
    ) %>%
    mutate(is_outlier = factor(is_outlier, levels = c(FALSE, TRUE))) %>%
    mutate(is_outlier_label = ifelse(is_outlier == FALSE, "Inlier", "Outlier")) %>%
    ungroup() %>%
    mutate(row_id = row_number())

  # Select top n lower and upper bound outliers
  top_lower_outliers <- data %>%
    filter(is_lower_outlier) %>%
    slice_min(FUSE_score, n = n_outliers, with_ties = FALSE)

  top_upper_outliers <- data %>%
    filter(is_upper_outlier) %>%
    slice_max(FUSE_score, n = n_outliers, with_ties = FALSE)

  # Combine top outliers
  top_outliers <- bind_rows(top_lower_outliers, top_upper_outliers)

  # Modify ClinVarLabelP for only top outliers
  data <- data %>%
    mutate(ClinVarLabelP = ifelse(row_id %in% unique(top_outliers$row_id),
      ClinVarLabelP, ""
    ))

  return(data)
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
    dplyr::arrange(variant_id)

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
