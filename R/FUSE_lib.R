library(DT)

df.funsum_all <- readRDS("~/github/ArchitectureOfCancer/data/funsum_maveDB_042423.rds")
df.NSFP_all <- read.table("~/github/ArchitectureOfCancer/data/dfNFSP_031022.txt", header = T, sep = "\t")
df.NSFP_all$aaref[which(df.NSFP_all$aaref == "X")] <- "*"
df.NSFP_all$aaalt[which(df.NSFP_all$aaalt == "X")] <- "*"


check_DMS_input <- function(df) {
  aa_list <- unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  msg <- c()
  if (ncol(df) < 5) {
    msg <- c(msg, "Error: DMS data has less than 5 columns!")
    return(list(df, msg))
  }

  col_names <- c("gene", "aapos", "aaref", "aaalt", "raw_score")
  if (!all(col_names %in% colnames(df))) {
    missing_col <- paste(col_names[!(col_names %in% colnames(df))], collapse = ",")
    msg <- c(msg, paste0("Error: DMS data missing ", missing_col, " columns!"))
    return(list(df, msg))
  }

  if ("Z" %in% df$aaref | "Z" %in% df$aaalt) {
    msg <- c(msg, "Warning: amino acid Z detected! Converted to stop codon '*'.")
    df$aaref[df$aaref == "Z"] <- "*"
    df$aaalt[df$aaalt == "Z"] <- "*"
  }

  if (!all(df$aaref %in% aa_list) | !all(df$aaalt %in% aa_list)) {
    msg <- c(msg, "Warning: unknown amino acid detected! Rows removed.")
    df <- df[df$aaref %in% aa_list, ]
    df <- df[df$aaalt %in% aa_list, ]
  }

  if (class(df$aapos) != "numeric") {
    msg <- c(msg, "Warning: non-numeric animo acid positions detected! Rows removed.")
    ind <- grep(pattern = "\\D", x = df$aapos)
    df <- df[-ind, ]
    df$aapos <- as.numeric(df$aapos)
  }

  if (class(df$raw_score) != "numeric") {
    msg <- c(msg, "Warning: non-numeric DMS score detected! Rows removed.")
    ind <- grep(pattern = "\\D", x = df$aapos)
    df <- df[-ind, ]
    df$aapos <- as.numeric(df$aapos)
  }

  return(list(df, msg))
}

check_direction <- function(df) {
  msg <- c()
  gene_id <- unique(df$gene)
  ind <- which((df$aaref != "P" & df$aaalt == "P") | (df$aaalt == "*"))
  if (length(ind) > 0) {
    if (mean(df$raw_score) > mean(df$raw_score[ind])) {
      df$raw_score <- -df$raw_score
      msg <- c(msg, paste0("Info: ", gene_id, " DMS data in opposite direction. Score inverted."))
    }
  } else {
    msg <- c(msg, paste0("Warning: ", gene_id, " DMS data direction cannot be determined!"))
  }

  return(list(df, msg))
}


# JS estimator
get_js <- function(m) {
  mbar <- mean(colMeans(m, na.rm = T), na.rm = T) # global mean
  mu0 <- colMeans(m, na.rm = T) # column means
  s2 <- var(as.vector(m), na.rm = T) / (nrow(m) * ncol(m) - length(which(is.na(as.vector(m)) == T))) # global variance
  cval <- 1 - (ncol(m) - 2) * s2 / sum((mu0 - mbar)^2, na.rm = T) # adjusted value
  js_est <- mbar + cval * (mu0 - mbar)
  return(js_est)
}

funsum_to_subTable <- function(df.funsum) {
  df.sub_tb <- c()
  for (i in 1:nrow(df.funsum)) {
    temp <- data.frame(aa_pair = paste0(rownames(df.funsum)[i], colnames(df.funsum)), score = as.vector(df.funsum[i, ]))
    df.sub_tb <- rbind(df.sub_tb, temp)
  }

  return(df.sub_tb)
}

de_noise <- function(df, pos_mean_method, df.funsum, include_LOF = T, show_func_class = F) {
  # filter out row without amino acid changes
  ind <- which(!is.na(df$aaref) & !is.na(df$aaalt))
  df <- df[ind, ]

  if (include_LOF) {
    aa_list <- unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  } else {
    aa_list <- unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
    ind <- which((df$aaref != "*") & (df$aaalt != "*"))
    df <- df[ind, ]
  }
  df.sub_tb <- funsum_to_subTable(df.funsum) # convert FUNSUM to tabular format

  if (is.null(df[["gene"]])) {
    df[["gene"]] <- "gene"
  }
  df$gene_aa_str <- paste0(df$gene, "---", df$aaref, df$aapos, df$aaalt)

  # collapse rows with the same amino acid substitutions
  df <- df %>%
    group_by(gene_aa_str) %>%
    reframe(
      gene = unique(gene), aapos = unique(aapos), aaref = unique(aaref), aaalt = unique(aaalt),
      raw_score = mean(raw_score), norm_raw_score = mean(norm_raw_score)
    )
  # summarise(
  #   gene = unique(gene), aapos = unique(aapos), aaref = unique(aaref), aaalt = unique(aaalt),
  #   raw_score = mean(raw_score), norm_raw_score = mean(norm_raw_score)
  # )
  
  ## calculate positional component
  df.pos_score <- df %>%
    group_by(gene, aapos) %>%
    reframe(aaref = unique(aaref), pos_mean = NA)
  #  summarise(aaref = unique(aaref), pos_mean = NA)
  df.pos_score$gene_aapos <- paste0(df.pos_score$gene, "---", df.pos_score$aapos)

  temp <- matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
  colnames(temp) <- aa_list
  if (pos_mean_method == "funsum") {
    temp2 <- matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
    colnames(temp2) <- aa_list
  }

  # pb = txtProgressBar(min = 1, max = nrow(df.pos_score), initial = 1, style = 3)
  for (i in 1:nrow(df.pos_score)) {
    ind <- which(df$gene == df.pos_score$gene[i] & df$aapos == df.pos_score$aapos[i])
    ind2 <- df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] <- df$norm_raw_score[ind[ind2]]

    # calculate pos_mean by funsum method
    if (pos_mean_method == "funsum") {
      ind <- which(!is.na(temp[i, ]))
      temp2[i, ind] <- temp[i, ind] - df.funsum[df.pos_score$aaref[i], aa_list[ind]]
    }

    # setTxtProgressBar(pb,i)
  }
  # close(pb)

  # calculate pos_mean by other methods
  if (pos_mean_method == "mean") {
    df.pos_score$pos_mean <- rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median") {
    df.pos_score$pos_mean <- apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js") {
    df.pos_score$pos_mean <- get_js(t(temp))
  } else if (pos_mean_method == "funsum") {
    df.pos_score$pos_mean <- get_js(t(temp2))
  }

  ## construct a new df with all possible substitutions
  df.out <- df.pos_score %>%
    dplyr::select(gene, aapos, aaref) %>%
    dplyr::slice(rep(1:n(), each = length(aa_list)))
  df.out$aaalt <- rep(aa_list, nrow(df.pos_score))

  # assign functional class
  if (show_func_class) {
    df.out$functional_class <- NA
    ind <- which(df.out$aaref != df.out$aaalt & df.out$aaalt != "*")
    df.out$functional_class[ind] <- "MIS"
    ind <- which(df.out$aaref != "*" & df.out$aaalt == "*")
    df.out$functional_class[ind] <- "LOF"
    ind <- which(df.out$aaref == df.out$aaalt)
    df.out$functional_class[ind] <- "SYN"
  }

  # assign norm_score, pos_score, sub_score to df.out
  df.out$gene_aa_str <- paste0(df.out$gene, "---", df.out$aaref, df.out$aapos, df.out$aaalt)
  df.out$gene_aapos <- paste0(df.out$gene, "---", df.out$aapos)
  df.out$aa_pair <- paste0(df.out$aaref, df.out$aaalt)
  ind <- match(df.out$gene_aa_str, table = df$gene_aa_str)
  df.out$raw_score <- df$raw_score[ind]
  df.out$norm_raw_score <- df$norm_raw_score[ind]
  ind <- match(df.out$gene_aapos, table = df.pos_score$gene_aapos)
  df.out$pos_score <- df.pos_score$pos_mean[ind]
  ind <- match(df.out$aa_pair, table = df.sub_tb$aa_pair)
  df.out$sub_score <- df.sub_tb$score[ind]
  df.out$final_score <- df.out$pos_score + df.out$sub_score
  df.out$final_score_lite <- df.out$final_score
  ind <- which(is.na(df.out$norm_raw_score))
  df.out$final_score_lite[ind] <- NA

  return(df.out)
}


get_FUSE_score_df <- function(filename, pos_mean_method = "js") {
  full_filename <- paste0(result_data_path, filename)
  df_all <- vroom(full_filename, delim = ",")
  temp <- check_DMS_input(df = df_all)
  df_all <- temp[[1]]
  msg <- temp[[2]]

  print(paste0("df_all nrow: ", nrow(df_all)))

  do_not_compute_fuse <- FALSE

  if (length(msg) > 0) {
    if (grepl(pattern = "error", x = msg, ignore.case = T)) {
      do_not_compute_fuse <- TRUE
    }
  }

  if (do_not_compute_fuse) {
    print("Failed data validity")
    return(DT())
  }

  gene_ids <- unique(df_all$gene)
  print(paste0("Info: number of genes: ", length(gene_ids)))

  ls.df <- list()
  df_all.out <- c()
  df.NSFP_all_out <- c()
  for (gene_id in gene_ids) {
    ## de-noise DMS data for each gene
    ind <- which(df_all$gene == gene_id)
    df <- df_all[ind, ]

    temp <- check_direction(df) # check and correct DMS data direction
    df <- temp[[1]]
    msg <- temp[[2]]
    print(msg)

    df.out <- de_noise(df, pos_mean_method = pos_mean_method, df.funsum = df.funsum_all)
    ls.df[[gene_id]] <- df.out
    df_all.out <- rbind(df_all.out, df.out)

    ## load df.NSFP for each gene
    ind <- which(df.NSFP_all$genename == gene_id)
    if (length(ind)) {
      df.NSFP <- df.NSFP_all[ind, ]
      df.NSFP$aa_str <- paste0(df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
      ind <- match(df.NSFP$aa_str, df.out$gene_aa_str)
      df.NSFP$norm_raw_score <- df.out$norm_raw_score[ind]
      df.NSFP$pos_score <- df.out$pos_score[ind]
      df.NSFP$sub_score <- df.out$sub_score[ind]
      df.NSFP$final_score <- df.out$final_score[ind]
      df.NSFP_all_out <- rbind(df.NSFP_all_out, df.NSFP)
    } else {
      print(paste0("Warning: no ClinVar record found for ", gene_id, "!"))
    }
  }

  ## render de-noised DMS data table
  ind <- which(df_all.out$aaref == df_all.out$aaalt)
  df_all.out <- df_all.out[-ind, c(
    "gene", "aapos", "aaref", "aaalt",
    "raw_score", "norm_raw_score", "pos_score", "sub_score", "final_score"
  )] %>%
    dplyr::filter(!(is.na(raw_score))) %>%
    dplyr::distinct()
  
  return(df_all.out)
}
