library(ggplot2)
library(cowplot)
library(rlang)
library(ggcorrplot)
library(RColorBrewer)
library(ggrepel)

getScatterPlot <- function(df,
                           title_txt,
                           xlabel,
                           ylabel,
                           x_col,
                           y_col,
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
  pt <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
      aes(color = ClinVarLabel),
      size = point_size,
      alpha = alpha
    ) +
    labs(x = xlabel, y = ylabel, title = title_txt) +
    theme_cowplot() +
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
    )

  if (annotate.point) {
    #pt <- pt + geom_text(aes(label = ClinVarLabelP), size = annotate_text_size) 
    pt <- pt + geom_text_repel(aes(label = ClinVarLabelP), size = annotate_text_size) 
    #      geom_text_repel(aes(label = ClinVarLabelP))
  } else {
    pt <- pt + geom_text_repel(aes(label = ClinVarLabel))
  }

  if (!is.null(xlimits) && !is.null(xbreaks)) {
    pt <- pt + scale_x_continuous(limits = xlimits, breaks = xbreaks)
  }
  if (!is.null(ylimits) && !is.null(ybreaks)) {
    pt <- pt + scale_y_continuous(limits = ylimits, breaks = ybreaks)
  }

  return(pt)
}

# getScatterPlot <- function(df,
#                            title_txt,
#                            xlabel,
#                            ylabel,
#                            x_col,
#                            y_col,
#                            alpha = 1,
#                            point_size = 1,
#                            point_color = "#00AFBB",
#                            title_font_size = 12,
#                            x_y_font_size = 12,
#                            annotate_text_size = 4,
#                            annotate.point = FALSE,
#                            xlimits = NULL,
#                            xbreaks = NULL,
#                            ylimits = NULL,
#                            ybreaks = NULL) {
#   pt <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
#     geom_point(
#       aes(color = ClinVarLabel),
#       size = point_size,
#       alpha = alpha
#     ) +
#     labs(x = xlabel, y = ylabel, title = title_txt) +
#     theme(
#       plot.title = element_text(
#         color = "black",
#         size = title_font_size,
#         face = "bold", hjust = 0.5
#       ),
#       axis.title.x = element_text(size = x_y_font_size, face = "bold"),
#       axis.title.y = element_text(size = x_y_font_size, face = "bold"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       axis.line = element_line(colour = "black")
#     ) +
#     theme_cowplot()
#
#   if (annotate.point) {
#     pt <- pt + geom_text(aes(label = ClinVarLabelP), size = annotate_text_size) +
#       geom_text_repel(aes(label = ClinVarLabelP))
#   } else {
#     pt <- pt +  geom_text_repel(aes(label = ClinVarLabel))
#   }
#
#   if (!is.null(xlimits) && !is.null(xbreaks)) {
#     pt <- pt + scale_x_continuous(limits = xlimits, breaks = xbreaks)
#   }
#   if (!is.null(ylimits) && !is.null(ybreaks)) {
#     pt <- pt + scale_y_continuous(limits = ylimits, breaks = ybreaks)
#   }
#
#   return(pt)
# }


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
