library(ggplot2)
library(cowplot)
library(rlang)
library(ggcorrplot)

get_scatter_plot <- function(df,
                             title_txt,
                             xlabel,
                             ylabel,
                             x_col,
                             y_col,
                             x_corr = 0,
                             y_corr = 10,
                             point_size = 1,
                             alpha = 1,
                             point_color = "#00AFBB",
                             font_size = 14,
                             annotate_text_size = 4,
                             xlimits = NULL,
                             xbreaks = NULL,
                             ylimits = NULL,
                             ybreaks = NULL) {
  # Compute correlation and p-value
  cor_test <- cor.test(
    df[[deparse(substitute(x_col))]],
    df[[deparse(substitute(y_col))]],
    method = "spearman"
  )
  cor_val <- cor_test$estimate
  p_val <- cor_test$p.value
  cor_text <- paste0(
    "Spearman Corr.: ",
    round(cor_val, 2), ",P-value: ", round(p_val, 4)
  )

  pt <- ggplot(df, aes(x = {{ x_col }}, y = {{ y_col }})) +
    geom_point(
      size = point_size,
      color = point_color,
      alpha = alpha
    ) +
    geom_smooth(method = lm) +
    labs(x = xlabel, y = ylabel, title = title_txt) +
    annotate("text",
      x = x_corr, y = y_corr, label = cor_text,
      hjust = .5, vjust = 1,
      size = annotate_text_size,
      colour = "darkcyan",
      fontface = "bold"
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
    theme_cowplot()

  if (!is.null(xlimits) && !is.null(xbreaks)) {
    pt <- pt + scale_x_continuous(limits = xlimits, breaks = xbreaks)
  }
  if (!is.null(ylimits) && !is.null(ybreaks)) {
    pt <- pt + scale_y_continuous(limits = ylimits, breaks = ybreaks)
  }

  pt
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
    theme_cowplot(12) +
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

  pt
}
