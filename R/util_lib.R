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
                                  ClinVarLabelP, ""))

  return(data)
}
