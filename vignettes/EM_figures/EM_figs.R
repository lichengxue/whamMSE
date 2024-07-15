# Load necessary libraries
library(ggplot2)
library(grid)
library(gridExtra)
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
    theme(plot.title = element_text(hjust = 0.4, size = 24, face = "bold")) +
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
      annotate("text", x = 6, y = y_position + 1.5, label = effects[1], color = effect_colors[1], size = 8, fontface = "bold") +
      annotate("segment", x = 5.8, xend = 6.2, y = y_position + 2, yend = y_position + 2, color = effect_colors[1], arrow = arrow(type = "closed")) +
      annotate("segment", x = 6.2, xend = 5.8, y = y_position + 1.0, yend = y_position + 1.0, color = effect_colors[1], arrow = arrow(type = "closed"))
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
          annotate("text", x = 6, y = y_position + 0.75, label = effect[eff_idx], color = effect_colors[[i]][eff_idx], size = 8, fontface = "bold") +
          annotate("segment", x = 5.8, xend = 6.2, y = y_position + 1.0, yend = y_position + 1.0, color = effect_colors[[i]][eff_idx], arrow = arrow(type = "closed")) +
          annotate("segment", x = 6.2, xend = 5.8, y = y_position + 0.5, yend = y_position + 0.5, color = effect_colors[[i]][eff_idx], arrow = arrow(type = "closed"))
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
vessel_positions_fig1 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7)) # Adjust these as needed
vessel_positions_fig2 <- list(c(4.5, 4.7, 3, 3.2), c(4.5, 4.7, 2, 2.2))
vessel_positions_fig3 <- list(c(4.5, 4.7, 1.5, 1.7)) # Only one vessel in fig.3
vessel_positions_fig4 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7))
vessel_positions_fig5 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7))
vessel_positions_fig6 <- list(c(4.5, 4.7, 1.5, 1.7), c(4.5, 4.7, 3.5, 3.7))

# Create figures
fig1 <- create_base_plot("Separate assessment models with NAA random effects",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("NAA RE"), c("NAA RE")),
                         list(c("royalblue"), c("green")),
                         vessel_positions = vessel_positions_fig1)

fig2 <- create_base_plot("One-area model with spatially-implicit fleets-as-areas\nand with NAA random effects",
                         c("Fleet 1", "Fleet 2"),
                         list(list(label = "Stock 1", color = "royalblue")),
                         list(c("NAA RE")),
                         list(c("royalblue")),
                         vessel_positions = vessel_positions_fig2,
                         overlap_stocks = TRUE,
                         split_box = TRUE)

fig3 <- create_base_plot("One-area model with catch aggregated\nand with NAA random effects",
                         c("Region 1"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Aggregated", color = "white")),
                         list(c("NAA RE")),
                         list(c("black")),
                         vessel_positions = vessel_positions_fig3)

# Create figure 4 with movement arrows
fig4 <- create_base_plot("Two-area model with fixed movement\nand with NAA random effects",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("NAA RE"), c("NAA RE")),
                         list(c("royalblue"), c("green")),
                         vessel_positions = vessel_positions_fig4) +
  geom_segment(aes(x = 5.5, y = 2.25, xend = 5.5, yend = 3.25),
               arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "black", size = 1) +
  geom_segment(aes(x = 4.5, y = 3.25, xend = 4.5, yend = 2.25),
               arrow = arrow(length = unit(0.2, "inches"), type = "closed"), color = "black", size = 1) +
  annotate("text", x = 5.6, y = 2.4, label = "Move", color = "black", size = 8, angle = 90, hjust = 0, fontface = "bold") +
  annotate("text", x = 4.4, y = 3, label = "Move", color = "black", size = 8, angle = 90, hjust = 1.1, fontface = "bold")


fig5 <- create_base_plot("Two-area model with no movement\nand with NAA random effects",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("NAA RE"), c("NAA RE")),
                         list(c("royalblue"), c("green")),
                         vessel_positions = vessel_positions_fig5)

fig6 <- create_base_plot("Two-area model with no movement\nand with Rec random effects",
                         c("Region 1", "Region 2"),
                         list(list(label = "Stock 1", color = "royalblue"), list(label = "Stock 2", color = "green")),
                         list(c("Rec RE"), c("Rec RE")),
                         list(c("royalblue"), c("green")),
                         vessel_positions = vessel_positions_fig6)

# Save plots
ggsave("Separate_Assessment_Models_with_NAA_Random_Effects.png", plot = fig1, width = 10, height = 6, bg = "transparent")
ggsave("One_Area_Model_with_Spatially_Implicit_Fleets_as_Areas_and_with_NAA_Random_Effects.png", plot = fig2, width = 10, height = 6, bg = "transparent")
ggsave("One_Area_Model_with_Catch_Aggregated_and_with_NAA_Random_Effects.png", plot = fig3, width = 10, height = 6, bg = "transparent")
ggsave("Two_Area_Model_with_Fixed_Movement_and_with_NAA_Random_Effects.png", plot = fig4, width = 10, height = 6, bg = "transparent")
ggsave("Two_Area_Model_with_No_Movement_and_with_NAA_Random_Effects.png", plot = fig5, width = 10, height = 6, bg = "transparent")
ggsave("Two_Area_Model_with_No_Movement_and_with_Rec_Random_Effects.png", plot = fig6, width = 10, height = 6, bg = "transparent")
