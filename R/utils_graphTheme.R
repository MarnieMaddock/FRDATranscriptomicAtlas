# utils_graphTheme.R
library(ggplot2)
theme_Marnie <- ggplot2::theme(axis.line.y = element_line(colour = "black", linewidth = 0.9),
                      axis.line.x = element_line(colour = "black", linewidth = 0.9),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_blank(),
                      axis.title.x = element_text(size = 18, margin = margin(5,0,0,0)),
                      axis.title.y = element_text(size = 18, margin = margin(0,10,0,0)),
                      axis.text = element_text(size = 18, colour = "black"),
                      axis.text.x = element_text(margin = margin(t=5), size=14),
                      axis.text.y = element_text(size=16),
                      plot.title = element_text(size = 30, hjust = 0), # legend.position = c(0.8, 0.8)
                      legend.position  = "right",
                      legend.key.size = unit(0.7, "cm"),
                      legend.text = element_text(size = 12),
                      legend.title = element_text(face = "bold", size = 14, hjust =0.5),
                      legend.key.width = unit(0.7,"cm"),
                      legend.key = element_rect(fill = NA, colour = NA),
                      strip.text = element_text(size = 16, face = "bold"),
                      strip.background = element_rect(colour = "black"),
                      panel.spacing = unit(0, "lines"),
                      plot.margin = margin(t=5, r=5, b=5, l=5, unit="pt") # Adjust these values accordingly

)
