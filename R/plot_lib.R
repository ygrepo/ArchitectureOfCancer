library(ggplot2)
library(cowplot)
library(rlang)
library(ggcorrplot)
library(RColorBrewer)
library(ggrepel)
library(ggthemes)

clinvar_label_color_mappings <- c(
  "None" = "#93C47D",
  "US" = "#6FA8DC",
  "LB/B" = "#CC00B1",
  "LP/P" = "#CC0000"
)

polyphen_label_color_mappings <- c(
  "benign" = "#6FA8DC",
  "possibly" = "#CC00B1",
  "probably" = "#CC0000"
)

beta_pval_cols <- c("up" = "red", "down" = "#26b3ff", "ns" = "grey")
beta_pval_sizes <- c("up" = 3, "down" = 3, "ns" = 2)
beta_pval_alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

theme_publication <- function(base_size = 12,
                              base_family = "Helvetica",
                              base_face = "bold",
                              legend_bottom = TRUE) {
  if (!is.null(legend_bottom)) {
    if (legend_bottom) {
      loc <- "bottom"
      orientation <- "horizontal"
    } else {
      loc <- "right"
      orientation <- "vertical"
    }

    return(ggthemes::theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
        plot.title = element_text(
          face = base_face,
          size = rel(1.2), hjust = 0.5
        ),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = base_face, size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = loc,
        legend.text = element_text(size = rel(1.2)),
        legend.direction = orientation,
        legend.key.size = unit(0.3, "cm"),
        legend.spacing = unit(0, "cm"),
        legend.title = element_text(face = "italic"),
        plot.margin = unit(c(5, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
        strip.text = element_text(face = base_face)
      ))
  }

  return(ggthemes::theme_foundation(base_size = base_size, base_family = base_family)
  + theme(
      plot.title = element_text(
        face = base_face,
        size = rel(1.2), hjust = 0.5
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = base_face, size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "none",
      plot.margin = unit(c(5, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = base_face)
    ))
}

getScatterPlot <- function(df,
                           title_txt,
                           xlabel,
                           ylabel,
                           x_col,
                           y_col,
                           size_col,
                           size_name,
                           legend_size_breaks,
                           legend_size_labels,
                           alpha = 1,
                           point_size = 1,
                           point_color = "#00AFBB",
                           title_font_size = 12,
                           x_y_font_size = 12,
                           annotate_text_size = 4,
                           annotate.point = FALSE,
                           xlimits = NULL,
                           xbreaks = NULL,
                           ylimits = NULL,
                           ybreaks = NULL) {
  size_col_sym <- rlang::sym(size_col)
  pt <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
      aes(
        color = ClinVar.Variant.Category,
        size = .data[[size_col]],
        alpha = alpha
      )
    ) +
    scale_size_continuous(
      name = size_name,
      range = c(1, 10),
      breaks = legend_size_breaks,
      labels = legend_size_labels
    ) +
    labs(x = xlabel, y = ylabel, title = title_txt) +
    guides(alpha = guide_none()) +
    theme(
      plot.title = element_text(
        color = "black",
        size = title_font_size,
        face = "bold", hjust = 0.5
      ),
      axis.title.x = element_text(size = x_y_font_size, face = "bold"),
      axis.title.y = element_text(size = x_y_font_size, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme_cowplot()

  if (annotate.point) {
    # pt <- pt + geom_text(aes(label = ClinVarLabelP), size = annotate_text_size)
    pt <- pt + geom_text_repel(aes(label = ClinVarLabelP),
      size = annotate_text_size,
      max.overlaps = 50
    )
    #      geom_text_repel(aes(label = ClinVarLabelP))
  } else {
    pt <- pt + geom_text_repel(aes(label = ClinVarLabel))
  }

  if (!is.null(xbreaks)) {
    pt <- pt + scale_x_continuous(breaks = xbreaks)
  }
  if (!is.null(xlimits)) {
    pt <- pt + xlim(xlimits)
  }

  if (!is.null(ybreaks)) {
    pt <- pt + scale_y_continuous(
      breaks = ybreaks
      #      expand = c(0.2, 0)
    )
  }
  if (!is.null(ylimits)) {
    pt <- pt + ylim(ylimits)
  }

  return(pt)
}


getCorrelationByCancerTypePlot <- function(df,
                                           lab_x_txt,
                                           expression_val,
                                           point_size = 3,
                                           point_alpha = 0.7,
                                           font_size = 14,
                                           legend_font_size = 12,
                                           xlimits = NULL,
                                           xbreaks = NULL) {
  pt <- ggplot(df, aes(x = Cor, y = Cancer, color = Study)) +
    geom_point(size = point_size, alpha = point_alpha) +
    labs(x = lab_x_txt, y = "Cancer") +
    labs(title = expression_val) +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(size = font_size, face = "bold"),
      axis.text.y = element_text(size = font_size, face = "bold"),
      axis.title.x = element_text(size = font_size, face = "bold"),
      axis.title.y = element_text(size = font_size, face = "bold"),
      plot.title = element_text(size = font_size, face = "bold"),
      legend.text = element_text(size = legend_font_size),
      legend.title = element_text(size = legend_font_size),
    )

  if (!is.null(xlimits) && !is.null(xbreaks)) {
    pt <- pt + scale_x_continuous(limits = xlimits, breaks = xbreaks)
  }

  return(pt)
}


get_violin_box_FUSE_score_by_pval_category <- function(df,
                                                       title_txt,
                                                       size_col,
                                                       xlabel,
                                                       ylabel,
                                                       title_font_size = 12,
                                                       x_y_font_size = 12,
                                                       annotate_text_size = 4,
                                                       alpha = 1,
                                                       annotate_flag = FALSE,
                                                       ybreaks = NULL) {
  pt <- ggplot(df, aes(x = pval_category, y = FUSE_score)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_point(
      aes(
        x = pval_category,
        y = FUSE_score,
        size = .data[[size_col]],
        color = ClinVarLabel,
        alpha = alpha,
        shape = is_outlier_label
      ),
      alpha = alpha,
      position = position_jitter(width = 0.2, height = 0)
    ) +
    scale_size_area() +
    scale_color_manual(values = clinvar_label_color_mappings) +
    scale_shape_manual(values = c(19, 17)) +
    labs(
      x = xlabel, y = ylabel, title = title_txt,
      color = "ClinVar Label",
      shape = "Inlier/Outlier"
    )

  if (annotate_flag) {
    pt <- pt + geom_text_repel(
      data = subset(df, is_outlier == TRUE),
      aes(
        x = pval_category,
        y = FUSE_score,
        label = ClinVarLabelP
      ),
      size = annotate_text_size,
      max.overlaps = 50
    )
  }
  pt <- pt +
    theme(
      plot.title = element_text(
        color = "black",
        size = title_font_size,
        face = "bold", hjust = 0.5
      ),
      axis.title.x = element_text(size = x_y_font_size, face = "bold"),
      axis.title.y = element_text(size = x_y_font_size, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme_cowplot()

  if (!is.null(ybreaks)) {
    pt <- pt + scale_y_continuous(
      breaks = ybreaks
    )
  }

  return(pt)
}

get_violin_box_polyphen_score_by_pval_category <- function(df,
                                                           p_val_df = NULL,
                                                           title_txt,
                                                           size_col,
                                                           var_label_col,
                                                           xlabel,
                                                           ylabel,
                                                           legend_bottom = FALSE,
                                                           title_font_size = 12,
                                                           x_y_font_size = 12,
                                                           annotate_text_size = 4,
                                                           alpha = 1,
                                                           annotate_flag = FALSE,
                                                           shape_manual = c(16, 18),
                                                           ybreaks = NULL) {
  pt <- ggplot(df, aes(x = pval_category, y = polyphen_score)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_point(
      aes(
        x = pval_category,
        y = polyphen_score,
        size = .data[[size_col]],
        color = polyphen_label,
        alpha = alpha,
        shape = is_outlier_label
      ),
      alpha = alpha,
      position = position_jitter(width = 0.2, height = 0)
    ) +
    scale_size_area() +
    scale_color_manual(values = polyphen_label_color_mappings) +
    scale_shape_manual(values = shape_manual) +
    labs(
      x = xlabel, y = ylabel, title = title_txt,
      color = "Polyphen Label",
      shape = "Inlier/Outlier"
    )

  if (annotate_flag) {
    pt <- pt + geom_text_repel(
      data = subset(df, is_outlier == TRUE),
      aes(
        x = pval_category,
        y = polyphen_score,
        label = .data[[var_label_col]]
      ),
      size = annotate_text_size,
      max.overlaps = 50
    )
  }

  if (!is.null(p_val_df)) {
    pt <- pt + ggprism::add_pvalue(p_val_df,
      xmin = "group1",
      xmax = "group2",
      label = "adjusted_p_value",
      y.position = "y.position",
      remove.bracket = TRUE
    )
    #             bracket.shorten = p_val_df$bracket_shorten)
  }

  pt <- pt + theme_publication(
    base_size = x_y_font_size,
    legend_bottom = legend_bottom
  ) +
    theme(
      plot.title = element_text(
        color = "black",
        size = title_font_size,
        face = "bold", hjust = 0.5
      ),
      axis.title.x = element_text(size = x_y_font_size, face = "bold"),
      axis.title.y = element_text(size = x_y_font_size, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme_cowplot()

  if (!is.null(ybreaks)) {
    pt <- pt + scale_y_continuous(
      breaks = ybreaks
    )
  }

  return(pt)
}


get_beta_pval_violin_plot <- function(df,
                                      title_txt,
                                      color_manual = beta_pval_cols,
                                      font_size = 14,
                                      legend_bottom = NULL) {
  sign_df <- df %>%
    dplyr::filter((beta_type == "down") |
      (beta_type == "up"))

  pt <- df %>%
    ggplot2::ggplot(aes(
      x = beta,
      y = LOG10ADJPVAL
    )) +
    geom_point(aes(colour = beta_type),
      alpha = 0.2,
      shape = 16,
      size = 1
    ) +
    geom_point(
      data = df %>%
        dplyr::filter(beta_type == "up"),
      shape = 21,
      size = 5,
      fill = "firebrick",
      colour = "black"
    ) +
    geom_point(
      data = df %>%
        dplyr::filter(beta_type == "down"),
      shape = 21,
      size = 5,
      fill = "steelblue",
      colour = "black"
    ) +
    labs(
      x = "Beta",
      y = "-log10(adj_p_value)",
      title = title_txt,
      color = "Beta"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c(-1, 2),
      linetype = "dashed"
    ) +
    ggrepel::geom_label_repel(
      data = sign_df,
      aes(label = ClinVarLabelP),
      force = 10,
      size = 4
      # max.overlaps = 50
    ) +
    theme_publication(
      base_size = font_size,
      legend_bottom = legend_bottom
    ) +
    theme(
      plot.title = element_text(
        color = "black",
        size = font_size,
        face = "bold", hjust = 0.5
      ),
      axis.title.x = element_text(size = font_size, face = "bold"),
      axis.title.y = element_text(size = font_size, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    scale_x_continuous(
      breaks = c(seq(-2, 6, 1)),
      limits = c(-2, 6)
    ) +
    scale_y_continuous(
      breaks = c(seq(0, 3, 0.5)),
      limits = c(0, 3)
    )

  return(pt)
}


get_pval_category_adj_pval_scatter_plot <- function(df,
                                                    title_txt,
                                                    font_size = 14,
                                                    legend_bottom = NULL,
                                                    max_overlaps_val = 50) {
  pt <- df %>%
    ggplot2::ggplot(aes(
      x = pval_cat,
      y = LOG10_adj_pval
    )) +
    geom_point(
      data = df,
      shape = 21,
      size = 5,
      fill = "steelblue",
      colour = "black"
    ) +
    labs(
      x = "Significance Test",
      y = "-log10(adj_p_value)",
      title = title_txt,
      color = "pval_cat"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    ggrepel::geom_label_repel(
      data = df,
      aes(label = gene),
      force = 10,
      size = 4,
      max.overlaps = max_overlaps_val
    ) +
    theme_publication(
      base_size = font_size,
      legend_bottom = legend_bottom
    ) +
    theme(
      plot.title = element_text(
        color = "black",
        size = font_size,
        face = "bold", hjust = 0.5
      ),
      axis.title.x = element_text(size = font_size, face = "bold"),
      axis.title.y = element_text(size = font_size, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    scale_y_continuous(
      breaks = c(seq(0, 1.5, 0.5)),
      limits = c(0, 1.5)
    )

  return(pt)
}


get_polyphen_score_genebass_pvalue <- function(df,
                                               title_txt,
                                               alpha = 1,
                                               title_font_size = 24,
                                               x_y_font_size = 24,
                                               legend_title_font_size = 14,
                                               legend_text_font_size = 8,
                                               stats_font_size = 6,
                                               regression_line_flag = FALSE) {
  pt <- df %>%
    ggplot2::ggplot(aes(
      x = polyphen_score,
      y = Log10PVal,
    )) +
    geom_point(
      aes(
        color = BETA,
        size = LOGAF
      ),
      alpha = alpha,
      shape = 16
    ) +
    labs(
      x = "Polyphen Score",
      y = "-log10(PValue)",
      size = "LOGAF",
      title = title_txt
    )

  if (regression_line_flag) {
    pt <- pt + ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      fullrange = TRUE,
      formula = y ~ x,
      linetype = "dashed",
      colour = "red"
    )
  }

  pt <- pt + geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
    cowplot::theme_cowplot() +
    theme(
      plot.title = element_text(
        color = "black",
        size = title_font_size,
        face = "bold",
        hjust = 0.5
      ),
      axis.title.x = element_text(size = x_y_font_size, face = "bold"),
      axis.title.y = element_text(size = x_y_font_size, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.text = element_text(size = legend_text_font_size),
      legend.title = element_text(face = "bold", size = legend_title_font_size)
    ) +
    scale_x_continuous(
      breaks = c(seq(0, 1.1, 0.2)),
      limits = c(0, 1.1)
    ) +
    scale_y_continuous(
      breaks = c(seq(0, 2.1, 0.2)),
      limits = c(0, 2.1)
    ) +
    scale_color_gradient2(
      low = "#bdbdbd",
      mid = "#e7e1ef",
      high = "#636363",
      name = "LOGAF"
    ) +
    scale_color_gradient2(
      low = "#c994c7",
      mid = "#e7e1ef",
      high = "#dd1c77",
      name = "BETA"
    )

  pt <- pt + ggpubr::stat_cor(
    method = "spearman",
    cor.coef.name = "R",
    size = stats_font_size,
    label.x = 0.6,
    label.y = 1,
    r.accuracy = 0.01,
    p.accuracy = 0.01
  )

  var_label_col <- "variantIdSign"
  pt <- pt + geom_text_repel(
    data = df,
    aes(
      x = polyphen_score,
      y = Log10PVal,
      label = .data[[var_label_col]]
    ),
    box.padding = .6,
    size = 6,
    max.overlaps = 10
  )

  return(pt)
}
