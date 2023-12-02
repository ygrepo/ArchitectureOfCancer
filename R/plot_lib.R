library(ggplot2)
library(cowplot)
library(rlang)
library(ggcorrplot)
library(RColorBrewer)
library(ggrepel)

clinvar_label_color_mappings <- c(
  "None" = "#93C47D",
  "US" = "#6FA8DC",
  "LB/B" = "#CC00B1",
  "LP/P" = "#CC0000"
)

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
    labs(x = xlabel, y = ylabel, title = title_txt, 
         color = "ClinVar Label",
         shape = "Inlier/Outlier")

  if (annotate_flag) {
    pt <- pt + geom_text_repel(
      data =subset(df, is_outlier == TRUE),
      aes(x = pval_category,
          y = FUSE_score,
          label = ClinVarLabelP),
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
