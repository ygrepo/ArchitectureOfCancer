library(dplyr)


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

get_genebass_polyphen_data_with_pval <- function(df,
                                                 description_val,
                                                 n_outliers = 5,
                                                 percent_IQR = 1.5) {
  df <- df %>%
    dplyr::filter(description == description_val) %>%
    dplyr::filter(!is.na(polyphen_score)) %>%
    dplyr::mutate(ProteinChange = str_extract(hgvsp, "(?<=:).+")) %>%
    dplyr::mutate(pval_category = case_when(
      genebass_pval <= 0.01 ~ "[0, 0.01]",
      (0.01 < genebass_pval) & (genebass_pval < 0.1) ~ "[0.01, 0.1[",
      (0.1 < genebass_pval) ~ "[0.1, 1[",
      TRUE ~ NA_character_
    )) %>%
    dplyr::mutate(pval_category = factor(pval_category)) %>%
    dplyr::mutate(LOG10AF = -log10(allele_frequency)) %>%
    dplyr::mutate(LOG10PVAL = -log10(genebass_pval)) %>%
    dplyr::mutate(ClinVarLabelP = gsub("p.", "", ProteinChange))


  df <- df %>%
    dplyr::filter(!is.na(pval_category)) %>%
    dplyr::group_by(pval_category) %>%
    dplyr::mutate(
      lower_bound = quantile(polyphen_score, 0.25, na.rm = TRUE)
      - percent_IQR * IQR(polyphen_score, na.rm = TRUE),
      upper_bound = quantile(polyphen_score, 0.75, na.rm = TRUE)
      + percent_IQR * IQR(polyphen_score, na.rm = TRUE),
      is_lower_outlier = ifelse(polyphen_score < lower_bound, TRUE, FALSE),
      is_upper_outlier = ifelse(polyphen_score > upper_bound, TRUE, FALSE),
      is_outlier = ifelse(is_lower_outlier | is_upper_outlier,
        TRUE, FALSE
      ),
    ) %>%
    dplyr::mutate(is_outlier = factor(is_outlier, levels = c(FALSE, TRUE))) %>%
    dplyr::mutate(is_outlier_label = ifelse(is_outlier == FALSE, "Inlier", "Outlier")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(row_id = row_number())

  # Select top n lower and upper bound outliers
  top_lower_outliers <- df %>%
    dplyr::filter(is_lower_outlier) %>%
    dplyr::slice_min(polyphen_score, n = n_outliers, with_ties = FALSE)

  top_upper_outliers <- df %>%
    dplyr::filter(is_upper_outlier) %>%
    dplyr::slice_max(polyphen_score, n = n_outliers, with_ties = FALSE)

  # Combine top outliers
  top_outliers <- dplyr::bind_rows(top_lower_outliers, top_upper_outliers)

  # Modify ClinVarLabelP for only top outliers
  df <- df %>%
    dplyr::mutate(ClinVarLabelP = ifelse(row_id %in% unique(top_outliers$row_id),
      ClinVarLabelP, ""
    ))
}


create_pval_category <- function(df,
                                 gene_val,
                                 description_val) {
  df <- df %>%
    dplyr::filter((description == description_val) & (gene == gene_val)) %>%
    dplyr::filter(!is.na(polyphen_score)) %>%
    dplyr::mutate(pval_category = case_when(
      Pvalue <= 0.01 ~ "[0, 0.01]",
      (0.01 < Pvalue) & (Pvalue <= 0.1) ~ "]0.01, 0.1]",
      (0.1 < Pvalue) ~ "]0.1, 1]",
      TRUE ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(pval_category)) %>%
    dplyr::mutate(pval_category = factor(pval_category))

  return(df)
}

process_polyphen_data_with_pval <- function(df,
                                            description_val,
                                            n_outliers = 5,
                                            percent_IQR = 1.5) {
  df <- df %>%
    dplyr::filter(description == description_val) %>%
    dplyr::filter(!is.na(polyphen_score)) %>%
    dplyr::mutate(pval_category = case_when(
      Pvalue <= 0.01 ~ "[0, 0.01]",
      (0.01 < Pvalue) & (Pvalue <= 0.1) ~ "]0.01, 0.1]",
      (0.1 < Pvalue) ~ "]0.1, 1]",
      TRUE ~ NA_character_
    )) %>%
    dplyr::mutate(pval_category = factor(pval_category)) %>%
    dplyr::mutate(LOG10AF = -log10(AF)) %>%
    dplyr::mutate(LOG10PVAL = -log10(Pvalue)) %>%
    dplyr::mutate(VarLabel = locus)


  df <- df %>%
    dplyr::filter(!is.na(pval_category)) %>%
    dplyr::group_by(pval_category) %>%
    dplyr::mutate(
      lower_bound = quantile(polyphen_score, 0.25, na.rm = TRUE)
      - percent_IQR * IQR(polyphen_score, na.rm = TRUE),
      upper_bound = quantile(polyphen_score, 0.75, na.rm = TRUE)
      + percent_IQR * IQR(polyphen_score, na.rm = TRUE),
      is_lower_outlier = ifelse(polyphen_score < lower_bound, TRUE, FALSE),
      is_upper_outlier = ifelse(polyphen_score > upper_bound, TRUE, FALSE),
      is_outlier = ifelse(is_lower_outlier | is_upper_outlier,
        TRUE, FALSE
      ),
    ) %>%
    dplyr::mutate(is_outlier = factor(is_outlier, levels = c(FALSE, TRUE))) %>%
    dplyr::mutate(is_outlier_label = ifelse(is_outlier == FALSE, "Inlier", "Outlier")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(row_id = row_number())

  # Select top n lower and upper bound outliers
  top_lower_outliers <- df %>%
    dplyr::group_by(pval_category) %>%
    dplyr::filter(is_lower_outlier) %>%
    dplyr::slice_min(polyphen_score, n = n_outliers, with_ties = FALSE)

  top_upper_outliers <- df %>%
    dplyr::group_by(pval_category) %>%
    dplyr::filter(is_upper_outlier) %>%
    dplyr::slice_max(polyphen_score, n = n_outliers, with_ties = FALSE)

  # Combine top outliers
  top_outliers <- dplyr::bind_rows(top_lower_outliers, top_upper_outliers)

  # Modify ClinVarLabelP for only top outliers
  df <- df %>%
    dplyr::mutate(VarLabel = ifelse(row_id %in% unique(top_outliers$row_id),
      VarLabel, ""
    ))

  return(df)
}

# Define a function to perform Wilcoxon rank-sum test and adjust p-values
get_wilcoxon_test_pairwise <- function(df1, df2, wilcox_method = "stats") {
  # print(paste0(group1," ",group2))
  if (wilcox_method == "stats") {
    wilcox_result <- stats::wilcox.test(df1, df2)
  }
  if (wilcox_method == "exactRankTest") {
    wilcox_result <- exactRankTests::wilcox.exact(df1, df2, exact = TRUE)
  }
  p_value <- wilcox_result$p.value
  # print(paste0("Pval:", p_value))
  return(p_value)
}


reorder_p_val_categories <- function(df) {
  df <- df %>%
    dplyr::mutate(temp = group1, group1 = group2, group2 = temp) %>%
    dplyr::select(-temp) %>%
    dplyr::slice(c(2, 3, 1)) %>%
    mutate(
      temp_group1 = case_when(
        group1 == "]0.1, 1]" & group2 == "]0.01, 0.1]" ~ group2,
        TRUE ~ group1
      ),
      temp_group2 = case_when(
        group1 == "]0.1, 1]" & group2 == "]0.01, 0.1]" ~ group1,
        TRUE ~ group2
      ),
      group1 = temp_group1,
      group2 = temp_group2
    ) %>%
    select(-temp_group1, -temp_group2)

  return(df)
}

perform_polyphen_score_signif_assoc_by_p_val_category <- function(df,
                                                                  cancer_type,
                                                                  gene_val,
                                                                  n_digits = 3) {
  df <- df %>%
    dplyr::filter((description == cancer_type) & (gene == gene_val))

  # Ensure pval_category is treated as a character vector
  df$pval_category <- as.character(df$pval_category)

  # Get unique levels of pval_category, excluding true NA values
  levels <- unique(df$pval_category)
  print(levels)
  # Initialize an empty dataframe to store the results
  pairwise_results <- tibble(
    cancer_type = character(),
    gene = character(),
    group1 = character(),
    group2 = character(),
    p_value = numeric()
  )

  if (length(levels) == 1) {
    print("Only one pval_category")
    pairwise_results <- tibble(
      cancer_type = cancer_type,
      gene = gene_val,
      group1 = NA,
      group2 = NA,
      p_value = NA
    )
    return(pairwise_results)
  }

  # Perform pairwise comparisons between levels of pval_category
  for (i in 1:(length(levels) - 1)) {
    for (j in (i + 1):length(levels)) {
      group1 <- levels[i]
      group2 <- levels[j]

      # print(paste0("i:", i, ",j:", j, "group1:", group1, ",group2:", group2))

      # Ensure group comparisons are only made for valid groups
      if ((length(group1) == 0) |
        is.na(group1) |
        (group1 == "") |
        (length(group2) == 0) |
        is.na(group2) |
        (group2 == "")) {
        print("Skip")
        next
      }
      if (group1 == group2) {
        print("Skip")
        next
      }
      df1 <- df$polyphen_score[df$pval_category == group1]
      df2 <- df$polyphen_score[df$pval_category == group2]

      p_value <- get_wilcoxon_test_pairwise(df1, df2)
      # print(paste0("Pval:", p_value))

      result_row <- tibble(
        cancer_type = cancer_type,
        gene = gene_val,
        group1 = group1,
        group2 = group2,
        p_value = p_value
      )

      pairwise_results <- rbind(pairwise_results, result_row)
    }
  }

  # Remove duplicated pairs, considering ordering
  pairwise_results <- pairwise_results %>%
    dplyr::distinct(group1, group2, .keep_all = TRUE)


  # Print and return the results
  # print(pairwise_results)

  return(reorder_p_val_categories(pairwise_results))
  #  return(reorder_p_val_categories(pairwise_results))
}
