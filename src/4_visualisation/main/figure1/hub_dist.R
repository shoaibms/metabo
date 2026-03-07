# -----------------------------------------------------------------------------
# Script: Hub Distribution and Rank Decay Analysis
#
# Description:
# This script generates two plots to visualize network hub properties from 
# metabolite data across different tissues (leaf, root) and genotypes (G1, G2).
#
# 1. Hub Distribution Plot: A violin and boxplot showing the distribution of
#    degree centrality for hubs.
# 2. Hub Rank Decay Plot: A line plot illustrating the decay in hub 
#    connectivity, ranked by degree.
#
# Input:
# - CSV files named 'leaf_g1_hub_metabolites.csv', 'leaf_g2_hub_metabolites.csv',
#   'root_g1_hub_metabolites.csv', 'root_g2_hub_metabolites.csv'.
#   Each file should contain at least a 'Degree' column.
#
# Output:
# - 'hub_distribution.pdf' and 'hub_distribution.png'
# - 'hub_rank_decay.pdf' and 'hub_rank_decay.png'
#
# -----------------------------------------------------------------------------

# --- 1. Load Libraries ---
library(ggplot2)
library(dplyr)

# --- 2. Configuration ---
# Hard-coded input directory as per user request
input_dir <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E"
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)

# Define color palettes for plots
palette_genotype <- c("G1" = "#2ecc71", "G2" = "#33d6d3")
palette_interaction <- c("Leaf.G1" = "#2ecc71", "Leaf.G2" = "#27ae60",
                         "Root.G1" = "#33d6d3", "Root.G2" = "#1e597d")

# --- 3. Data Loading and Preparation ---
# Load individual datasets
files <- list(
  leaf_g1 = read.csv(file.path(input_dir, "leaf_g1_hub_metabolites.csv")),
  leaf_g2 = read.csv(file.path(input_dir, "leaf_g2_hub_metabolites.csv")), 
  root_g1 = read.csv(file.path(input_dir, "root_g1_hub_metabolites.csv")),
  root_g2 = read.csv(file.path(input_dir, "root_g2_hub_metabolites.csv"))
)

# Create a single, combined dataframe with tissue and genotype labels
combined_df <- bind_rows(
  mutate(files$leaf_g1, Tissue = "Leaf", Genotype = "G1"),
  mutate(files$leaf_g2, Tissue = "Leaf", Genotype = "G2"),
  mutate(files$root_g1, Tissue = "Root", Genotype = "G1"),
  mutate(files$root_g2, Tissue = "Root", Genotype = "G2")
)

# --- 4. Plotting ---
# Define a custom theme for consistent plot appearance
theme_custom <- theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Plot 1: Hub Distribution
hub_dist_plot <- ggplot(combined_df, aes(x = Tissue, y = Degree, fill = Genotype)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.7, position = position_dodge(0.9)) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.3, position = position_jitterdodge(dodge.width=0.9)) +
  stat_summary(fun.data = "mean_sdl", geom = "pointrange", position = position_dodge(0.9)) +
  scale_fill_manual(values = palette_genotype) +
  theme_custom +
  labs(
    title = "Hub Molecular Features Distribution",
    subtitle = "Degree centrality across tissues and genotypes",
    y = "Degree Centrality",
    x = "Tissue"
  )

# Plot 2: Hub Rank Decay
# Prepare data by sorting by degree and adding a rank column for each group
rank_decay_df <- combined_df %>%
  arrange(desc(Degree)) %>%
  group_by(Tissue, Genotype) %>%
  mutate(Rank = row_number())

rank_decay_plot <- ggplot(rank_decay_df, aes(x = Rank, y = Degree, color = interaction(Tissue, Genotype))) +
  geom_line(size = 1) +
  geom_smooth(se = TRUE, alpha = 0.2, method = "loess", formula = y ~ x) +
  scale_color_manual(values = palette_interaction) +
  theme_custom +
  labs(
    title = "Hub Connectivity Decay",
    subtitle = "Degree distribution by rank order",
    x = "Hub Rank (ordered by degree)",
    y = "Degree Centrality",
    color = "Tissue-Genotype"
  )

# --- 5. Save Plots ---
# Function to save plots in multiple formats
save_plot <- function(plot, name, output_dir, width, height) {
  formats <- c("pdf", "png")
  for (format in formats) {
    filename <- file.path(output_dir, paste0(name, ".", format))
    ggsave(
      filename, 
      plot = plot, 
      width = width, 
      height = height, 
      dpi = 300
    )
  }
}

# Save the generated plots
save_plot(hub_dist_plot, "hub_distribution", output_dir, 6, 4)
save_plot(rank_decay_plot, "hub_rank_decay", output_dir, 6, 4)
