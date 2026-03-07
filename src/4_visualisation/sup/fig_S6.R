
library(ggplot2)
library(dplyr)
library(tidyr)

# Set output directory
output_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read data
data <- read.csv("C:/Users/ms/Desktop/r/chem_data/metabo2/result/combine_result/temporal_coherence.csv")

# Prepare data - clean pathway names and convert coherence to factor
plot_data <- data %>%
  mutate(
    Pathway_Clean = gsub(" metabolism", "", 
                         gsub("Biosynthesis of ", "",
                              gsub("Carbon fixation by ", "", Pathway_Name))),
    Tissue_Genotype = factor(paste(Tissue.type, Genotype),
                             levels = c("L G1", "L G2", "R G1", "R G2")),
    # Convert coherence to factor with clear labels
    coherence_factor = factor(temporal_coherence,
                              levels = c(-1, -0.5, 0.5, 1),
                              labels = c("Strong Negative (-1.0)",
                                         "Moderate Negative (-0.5)",
                                         "Moderate Positive (0.5)",
                                         "Strong Positive (1.0)"))
  )

# Create main heatmap with discrete colors
p <- ggplot(plot_data, 
            aes(x = Tissue_Genotype, 
                y = reorder(Pathway_Clean, temporal_coherence), 
                fill = coherence_factor)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Strong Negative (-1.0)" = "#1b8778",
               "Moderate Negative (-0.5)" = "#60ccbd",
               "Moderate Positive (0.5)" = "#e2f29b",
               "Strong Positive (1.0)" = "#abbd5e"),
    name = "Temporal\nCoherence"
  ) +
  labs(
    x = "Tissue-Genotype",
    y = "Metabolic Pathway",
    title = "Metabolic Pathway Temporal Coherence"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

# Save plots
ggsave(
  filename = file.path(output_dir, "temporal_coherence_heatmap.pdf"),
  plot = p,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "temporal_coherence_heatmap.png"),
  plot = p,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

print("Plots saved in result/section3 directory")

