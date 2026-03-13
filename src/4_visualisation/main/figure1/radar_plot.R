#-------------------------------------------------------------------------------
#
#                           Radar Plot of Network Metrics
#
#-------------------------------------------------------------------------------
#
# Description:
#   This script generates a radar plot to visualize and compare network metrics
#   across different tissue types and genotypes. The plot is saved as both PDF
#   and PNG files.
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  1. Load required libraries
#-------------------------------------------------------------------------------
# It is a good practice to load all necessary packages at the start of a script.
library(ggplot2)
library(dplyr)
library(ggradar)
library(scales) # ggradar may use this under the hood.

#-------------------------------------------------------------------------------
#  2. Set file paths
#-------------------------------------------------------------------------------
# NOTE: For publication, it is recommended to use relative paths or command-line
# arguments for input/output files to ensure reproducibility.
# File paths
input_path <- "C:/Users/USER/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/network_metrics_summary.csv"
output_dir <- "C:/Users/USER/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/plots"

# Create the output directory if it does not exist.
# `showWarnings = FALSE` prevents a warning if the directory already exists.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
#  3. Read and preprocess data
#-------------------------------------------------------------------------------
# Read the network metrics data from the specified CSV file.
network_metrics <- read.csv(input_path)

# Prepare the data for the radar plot.
# - A 'group' column is created by combining 'Tissue.type' and 'Genotype'.
# - Relevant metric columns are selected.
radar_data <- network_metrics %>%
  mutate(group = paste(Tissue.type, Genotype)) %>%
  select(
    group,
    density,
    transitivity,
    modularity,
    temporal_coherence,
    path_length,
    directionality
  )

# Scale the data for all metrics to a range of 0 to 1.
# This normalization is necessary for plotting different metrics on the same radar chart axes.
# `across()` applies the scaling function to all columns except 'group'.
radar_data_scaled <- radar_data %>%
  mutate(across(-group, ~ rescale(.))) # Using scales::rescale for clarity

# Rename columns to be more descriptive for the plot labels.
radar_data_formatted <- radar_data_scaled %>%
  rename(
    "Density" = density,
    "Transitivity" = transitivity,
    "Modularity" = modularity,
    "Temporal Coherence" = temporal_coherence,
    "Path Length" = path_length,
    "Directionality" = directionality
  )

#-------------------------------------------------------------------------------
#  4. Generate and customize the radar plot
#-------------------------------------------------------------------------------
# Create the radar plot using the `ggradar` function.
radar_plot <- ggradar(
  radar_data_formatted,
  grid.min = 0,
  grid.mid = 0.5,
  grid.max = 1,
  values.radar = c("0", "0.5", "1"),
  group.point.size = 3,
  group.line.width = 1,
  # Using a standard color palette is good practice.
  # These colors are chosen for good visibility.
  group.colours = c("#2ecc71", "#27ae60", "#33d6d3", "#1e597d"),
  background.circle.colour = "white",
  gridline.min.colour = "grey90",
  gridline.mid.colour = "grey85",
  gridline.max.colour = "grey80",
  legend.position = "top",
  legend.title = "Tissue",
  plot.title = "Network Metrics Comparison",
  fill = TRUE,
  fill.alpha = 0.25,
  axis.label.size = 7,
  grid.label.size = 8,
  plot.extent.x.sf = 1.2, # Adjust plot extension to prevent label cutoff
  plot.extent.y.sf = 1.2
) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), # Center title
    legend.position = "top",
    legend.justification = "center",
    plot.margin = margin(5, 10, 5, 10),
    legend.spacing.x = unit(0.3, 'cm'),
    legend.margin = margin(0, 0, 10, 0)
  )

#-------------------------------------------------------------------------------
#  5. Save the plot to files
#-------------------------------------------------------------------------------
# Save the plot in both PDF (vector) and PNG (raster) formats for publication.
ggsave(
  file.path(output_dir, "network_metrics_radar.pdf"),
  plot = radar_plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(output_dir, "network_metrics_radar.png"),
  plot = radar_plot,
  width = 7,
  height = 6,
  dpi = 300
)
