# Load necessary libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(png)

# Function to create an ellipse
create_ellipse <- function(cx, cy, rx, ry, npoints = 100) {
  theta <- seq(0, 2 * pi, length.out = npoints)
  x <- cx + rx * cos(theta)
  y <- cy + ry * sin(theta)
  data.frame(x, y)
}

# Function to load an image
load_image <- function(image_path) {
  rasterGrob(readPNG(image_path), width = unit(2, "npc"), height = unit(3, "npc"))
}

# Function to create a base plot with fishing vessel
create_base_plot <- function(title, regions, stocks, effects, effect_colors, vessel_positions, overlap_stocks = FALSE, split_box = FALSE) {
  p <- ggplot() + 
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.4, size = 18, face = "bold")) +
    ggtitle(title)
  
  y_position <- 1
  
  if (split_box) {
    p <- p + 
      annotate("rect", xmin = 4, xmax = 6, ymin = y_position, ymax = y_position + 3, alpha = 1, color = "black", linetype = "dashed", fill = "white") +
      annotate("segment", x = 4, xend = 6, y = y_position + 1.5, yend = y_position + 1.5, color = "black", linetype = "dashed") +
      annotate("text", x = 5, y = y_position + 2.25, label = regions[1], hjust = 0.5, size = 8, fontface = "bold") +
      annotate("text", x = 5, y = y_position + 0.75, label = regions[2], hjust = 0.5, size = 8, fontface = "bold") +
      geom_polygon(data = create_ellipse(5, y_position + 1.5, 0.5, 1), aes(x, y), fill = stocks[[1]]$color, color = "black", alpha = 0.5) +
      annotate("text", x = 5, y = y_position + 1.5, label = stocks[[1]]$label, hjust = 0.5, size = 8, fontface = "bold") +
      annotate("text", x = 6, y = y_position + 1.5, label = effects[1], color = "black", size = 8, fontface = "bold") +
      annotate("segment", x = 5.8, xend = 6.2, y = y_position + 2, yend = y_position + 2, color = "black", arrow = arrow(type = "closed")) +
      annotate("segment", x = 6.2, xend = 5.8, y = y_position + 1.0, yend = y_position + 1.0, color = "black", arrow = arrow(type = "closed"))
  } else {
    for (i in 1:length(regions)) {
      p <- p + 
        annotate("rect", xmin = 4, xmax = 6, ymin = y_position, ymax = y_position + 1.5, alpha = 1, color = "black", linetype = "dashed", fill = "white") +
        annotate("text", x = 4.1, y = y_position + 0.75, label = regions[i], hjust = 0, size = 8, fontface = "bold")
      
      if (!overlap_stocks) {
        stock <- stocks[[i]]
        ellipse_data <- create_ellipse(5, y_position + 0.75, 0.5, 0.5)
        p <- p + 
          geom_polygon(data = ellipse_data, aes(x, y), fill = stock$color, color = "black", alpha = 0.5) +
          annotate("text", x = 5, y = y_position + 0.75, label = stock$label, hjust = 0.5, size = 8, fontface = "bold")
      }
      
      effect <- effects[[i]]
      for (eff_idx in 1:length(effect)) {
        p <- p + 
          annotate("text", x = 6, y = y_position + 0.75, label = effect[eff_idx], color = "black", size = 8, fontface = "bold") +
          annotate("segment", x = 5.8, xend = 6.2, y = y_position + 1.0, yend = y_position + 1.0, color = "black", arrow = arrow(type = "closed")) +
          annotate("segment", x = 6.2, xend = 5.8, y = y_position + 0.5, yend = y_position + 0.5, color = "black", arrow = arrow(type = "closed"))
      }
      
      y_position <- y_position + 2
    }
  }
  
  # Add fishing vessel image
  for (pos in vessel_positions) {
    p <- p + annotation_custom(load_image("vessel.png"), xmin = pos[1], xmax = pos[2], ymin = pos[3], ymax = pos[4])
  }
  
  return(p)
}

# Define vessel positions for each figure
vessel_positions_fig1 <- list(c(4.5, 4.7, 1.5, 1.7)) # Only one vessel in fig.1
vessel_positions_fig2 <- list(c(4.5, 4.7, 3, 3.2), c(4.5, 4.7, 2, 2.2))
vessel_positions_fig3 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7)) # 
vessel_positions_fig4 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7))
vessel_positions_fig5 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7))
vessel_positions_fig6 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7))

# Create figures with updated titles and natal homing arrows
fig1 <- create_base_plot("EM1: Panmictic (catch aggregated)\nNo movement, NAA random effects",
                         c("Region 1"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Aggregated", color = "white")),
                         list(c("NAA RE")),
                         list(c("black")),
                         vessel_positions = vessel_positions_fig1)

fig2 <- create_base_plot("EM2: Spatially implicit (fleets-as-areas)\nNo movement, NAA random effects",
                         c("Fleet 1", "Fleet 2"),
                         list(list(label = "Stock 1", color = "royalblue")),
                         list(c("NAA RE")),
                         list(c("black")),
                         vessel_positions = vessel_positions_fig2,
                         overlap_stocks = TRUE,
                         split_box = TRUE)

fig3 <- create_base_plot("EM3: Separate panmictic\nNo movement, NAA random effects",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("NAA RE"), c("NAA RE")),
                         list(c("black"), c("black")),
                         vessel_positions = vessel_positions_fig3)

fig4 <- create_base_plot("EM4: Spatially disaggregated\nNo movement, NAA random effects, Global SPR-based Brps",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("NAA RE"), c("NAA RE")),
                         list(c("black"), c("black")),
                         vessel_positions = vessel_positions_fig4)

fig5 <- create_base_plot("EM5: Spatially explicit\nFixed movement, NAA and movement random effects, Global SPR-based Brps",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("NAA RE"), c("NAA RE")),
                         list(c("black"), c("black")),
                         vessel_positions = vessel_positions_fig5) +
  geom_segment(aes(x = 5.5, y = 2.25, xend = 5.5, yend = 3.25),
               arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "black", size = 1) +
  geom_segment(aes(x = 4.5, y = 3.25, xend = 4.5, yend = 2.25),
               arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "black", size = 1) +
  # Adding natal homing arrows
  geom_curve(aes(x = 4.7, y = 2.2, xend = 4.7, yend = 3.3),
             arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "blue", size = 1, curvature = -0.3) +
  geom_curve(aes(x = 5.3, y = 3.3, xend = 5.3, yend = 2.2),
             arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "blue", size = 1, curvature = -0.3) +
  annotate("text", x = 5.6, y = 2.4, label = "Move", color = "black", size = 8, angle = 90, hjust = 0, fontface = "bold") +
  annotate("text", x = 4.4, y = 3, label = "Move", color = "black", size = 8, angle = 90, hjust = 1.1, fontface = "bold") +
  annotate("text", x = 5, y = 2.8, label = "Natal Homing", color = "blue", size = 6, hjust = 0.5, fontface = "bold") 

fig6 <- create_base_plot("EM6: Spatially explicit\nEstimated movement, NAA and movement random effects, Global SPR-based Brps",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("NAA RE"), c("NAA RE")),
                         list(c("black"), c("black")),
                         vessel_positions = vessel_positions_fig6) +
  geom_segment(aes(x = 5.5, y = 2.25, xend = 5.5, yend = 3.25),
               arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "black", size = 1, linetype = "dashed") +
  geom_segment(aes(x = 4.5, y = 3.25, xend = 4.5, yend = 2.25),
               arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "black", size = 1, linetype = "dashed") +
  # Adding natal homing arrows
  geom_curve(aes(x = 4.7, y = 2.2, xend = 4.7, yend = 3.3),
             arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "blue", size = 1, curvature = -0.3) +
  geom_curve(aes(x = 5.3, y = 3.3, xend = 5.3, yend = 2.2),
             arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "blue", size = 1, curvature = -0.3) +
  annotate("text", x = 5.6, y = 2.4, label = "Move", color = "black", size = 8, angle = 90, hjust = 0, fontface = "bold") +
  annotate("text", x = 4.4, y = 3, label = "Move", color = "black", size = 8, angle = 90, hjust = 1.1, fontface = "bold") +
  annotate("text", x = 5, y = 2.8, label = "Natal Homing", color = "blue", size = 6, hjust = 0.5, fontface = "bold") 

# Save individual plots
ggsave("EM1_Panmictic_Catch_Aggregated.png", plot = fig1, width = 10, height = 6, bg = "transparent")
ggsave("EM2_Spatially_Implicit_Fleets_as_Areas.png", plot = fig2, width = 10, height = 6, bg = "transparent")
ggsave("EM3_Separate_Panmictic.png", plot = fig3, width = 10, height = 6, bg = "transparent")
ggsave("EM4_Spatially_Disaggregated_No_Movement.png", plot = fig4, width = 10, height = 6, bg = "transparent")
ggsave("EM5_Spatially_Explicit_Fixed_Movement.png", plot = fig5, width = 10, height = 6, bg = "transparent")
ggsave("EM6_Spatially_Explicit_Estimated_Movement.png", plot = fig6, width = 10, height = 6, bg = "transparent")

# Add a thicker border to each figure using theme()
border_size <- 2
fig1 <- fig1 + theme(plot.background = element_rect(color = "black", size = border_size))
fig2 <- fig2 + theme(plot.background = element_rect(color = "black", size = border_size))
fig3 <- fig3 + theme(plot.background = element_rect(color = "black", size = border_size))
fig4 <- fig4 + theme(plot.background = element_rect(color = "black", size = border_size))
fig5 <- fig5 + theme(plot.background = element_rect(color = "black", size = border_size))
fig6 <- fig6 + theme(plot.background = element_rect(color = "black", size = border_size))

# Combine all figures into one plot arranged in 3 rows x 2 columns
combined_plot <- plot_grid(fig1, fig2, fig3, fig4, fig5, fig6, ncol = 2, nrow = 3)

# Save the combined plot as a single PNG
ggsave("Combined_Figures.png", plot = combined_plot, width = 20, height = 18, bg = "transparent")
