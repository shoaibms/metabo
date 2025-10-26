#-------------------------------------------------------------------------------#
#                                                                               #
#         Hub Distribution and Hub Rank Decay Plotting Script                   #
#                                                                               #
#-------------------------------------------------------------------------------#

# Description:
# This script generates two plots to visualize network hub properties from 
# metabolomics data:
# 1. Hub Distribution Plot: A violin and boxplot showing the distribution of 
#    degree centrality for molecular features across different tissues and genotypes.
# 2. Hub Rank Decay Plot: A line plot showing the decay of degree centrality 
#    as a function of hub rank, illustrating the connectivity distribution.
#
# The script reads hub metabolite data, prepares it for plotting, creates the
# visualizations using ggplot2, and saves them as PDF and PNG files.

#--- 1. Load Libraries ---#
library(ggplot2)
library(dplyr)

#--- 2. Configuration ---#

# NOTE: The input directory is hardcoded as requested. For better reproducibility
# and sharing, consider using relative paths or a project management tool 
# like the 'here' package in the future.
input_dir <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E"
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)

#--- 3. Load and Prepare Data ---#

# List hub metabolite files
files_to_read <- list(
  leaf_g1 = file.path(input_dir, "leaf_g1_hub_metabolites.csv"),
  leaf_g2 = file.path(input_dir, "leaf_g2_hub_metabolites.csv"), 
  root_g1 = file.path(input_dir, "root_g1_hub_metabolites.csv"),
  root_g2 = file.path(input_dir, "root_g2_hub_metabolites.csv")
)

# Check for file existence before reading
files_exist <- sapply(files_to_read, file.exists)
if (!all(files_exist)) {
  stop("Error: Not all input files were found. Missing: ", 
       paste(names(files_to_read)[!files_exist], collapse = ", "))
}

data_list <- lapply(files_to_read, read.csv)

# Combine datasets and add identifying columns
combined_df <- bind_rows(
  mutate(data_list$leaf_g1, Tissue = "Leaf", Genotype = "G1"),
  mutate(data_list$leaf_g2, Tissue = "Leaf", Genotype = "G2"),
  mutate(data_list$root_g1, Tissue = "Root", Genotype = "G1"),
  mutate(data_list$root_g2, Tissue = "Root", Genotype = "G2")
)

#--- 4. Define Plotting Theme ---#

# A consistent theme for publication-quality plots
publication_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

#--- 5. Create Plots ---#

# 5.1. Hub Distribution Plot
hub_dist_plot <- ggplot(combined_df, aes(x = Tissue, y = Degree, fill = Genotype)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.7, position = position_dodge(0.9)) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.3, position = position_jitterdodge(dodge.width=0.9)) +
  scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
  stat_summary(fun.data = "mean_sdl", geom = "pointrange", 
               position = position_dodge(0.9), color = "black") +
  publication_theme +
  labs(
    title = "Hub Molecular Features Distribution",
    subtitle = "Degree centrality across tissues and genotypes",
    y = "Degree Centrality", 
    x = "Tissue"
  )

# 5.2. Hub Rank Decay Plot

# Prepare data for rank decay plot: rank hubs within each group
ranked_df <- combined_df %>%
  group_by(Tissue, Genotype) %>%
  arrange(desc(Degree)) %>%
  mutate(Hub_Rank = row_number()) %>%
  ungroup()

rank_decay_plot <- ggplot(ranked_df, aes(x = Hub_Rank, y = Degree, color = interaction(Tissue, Genotype))) +
  geom_line(size = 1) +
  geom_smooth(se = TRUE, alpha = 0.2, method = "loess", formula = y ~ x) +
  scale_color_manual(
    name = "Tissue-Genotype",
    values = c(
      "Leaf.G1" = "#2ecc71", 
      "Leaf.G2" = "#27ae60",
      "Root.G1" = "#33d6d3", 
      "Root.G2" = "#1e597d"
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  publication_theme +
  labs(
    title = "Hub Connectivity Decay",
    subtitle = "Degree distribution by rank order",
    x = "Hub Rank (ordered by degree)", 
    y = "Degree Centrality"
  )

#--- 6. Save Plots ---#

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
save_plot(hub_dist_plot, "hub_distribution", output_dir, 8, 6)
save_plot(rank_decay_plot, "hub_rank_decay", output_dir, 8, 6)

#--- Script End ---#
