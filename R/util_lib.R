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
    mutate(ProteinChange = str_extract(ProteinChange, "(?<=:).+"))

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
