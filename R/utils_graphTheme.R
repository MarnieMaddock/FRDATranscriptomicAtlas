#' Package theme
#' @importFrom grid unit
#' @importFrom ggplot2 margin
#' @noRd
theme_Marnie <- function(base_size = 14) {
  ggplot2::theme(
    axis.line.y      = ggplot2::element_line(colour = "black", linewidth = 0.9),
    axis.line.x      = ggplot2::element_line(colour = "black", linewidth = 0.9),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white"),
    panel.border     = ggplot2::element_blank(),

    axis.title.x = ggplot2::element_text(
      size = base_size + 4,
      margin = ggplot2::margin(5, 0, 0, 0)
    ),
    axis.title.y = ggplot2::element_text(
      size = base_size + 4,
      margin = ggplot2::margin(0, 10, 0, 0)
    ),

    axis.text   = ggplot2::element_text(size = base_size, colour = "black"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
    axis.text.y = ggplot2::element_text(),

    plot.title = ggplot2::element_text(size = base_size + 16, hjust = 0),

    legend.position  = "right",
    legend.key.size  = grid::unit(0.7, "cm"),
    legend.text      = ggplot2::element_text(size = base_size - 2),
    legend.title     = ggplot2::element_text(face = "bold", size = base_size, hjust = 0.5),
    legend.key.width = grid::unit(0.7, "cm"),
    legend.key       = ggplot2::element_rect(fill = NA, colour = NA),

    strip.text       = ggplot2::element_text(size = base_size + 2, face = "bold"),
    strip.background = ggplot2::element_rect(colour = "black"),

    panel.spacing = grid::unit(0, "lines"),
    plot.margin   = ggplot2::margin(5, 5, 5, 5, unit = "pt")
  )
}

