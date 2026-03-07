################################################################################
#
# Module-Level Organisation Analysis and Visualization
#
# Author: [Your Name/Organization]
# Date: [Date]
#
# Description:
# This script analyzes the module structure of biological networks from different
# tissues (leaf and root) and genotypes. It uses the Louvain community detection
# algorithm to identify modules, calculates module sizes, and generates a plot
# to visualize the module organization across tissues and genotypes.
#
# Input:
#   - A '.rds' file containing network results, expected to have 'leaf_network'
#     and 'root_network' objects.
#     (Path: C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000/network_results.rds)
#
# Output:
#   - A summary statistics file ('module_summary_statistics.txt').
#   - Publication-quality plots of module organization in PDF, TIFF, and PNG
#     formats ('module_organization.pdf', 'module_organization.tiff',
#     'module_organization.png').
#
# Dependencies:
#   - igraph, ggplot2, tidyverse, gridExtra
#   - Ensure these packages are installed: install.packages(c("igraph", "ggplot2", "tidyverse", "gridExtra"))
#
################################################################################

# Load required libraries
library(igraph)
library(ggplot2)
library(tidyverse)

# --- Configuration ---

# Hard-coded base directory for input data as requested.
base_dir <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"
output_dir <- file.path(base_dir, "module_analysis")

# Create the output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Constants
EDGE_PROBABILITY_THRESHOLD <- 0.5

# --- Function Definitions ---

#' Analyze Network Modules
#'
#' Detects communities in a network using the Louvain algorithm and returns
#' module assignments for each node.
#'
#' @param network_data A list containing an 'edge_prob' matrix.
#' @param tissue_type A string indicating the tissue type (e.g., "Leaf").
#' @param genotype A string indicating the genotype (e.g., "G1").
#' @return A data frame with node, module ID, module size, tissue, and genotype.
analyze_modules <- function(network_data, tissue_type, genotype) {
  adj_matrix <- network_data$edge_prob > EDGE_PROBABILITY_THRESHOLD
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  
  louvain_comm <- cluster_louvain(g)
  module_sizes <- sizes(louvain_comm)
  
  module_assignments <- data.frame(
    node = V(g)$name,
    module_id = membership(louvain_comm),
    module_size = module_sizes[membership(louvain_comm)],
    tissue = tissue_type,
    genotype = genotype
  )
  
  return(module_assignments)
}

#' Save Publication-Quality Plots
#'
#' Saves a ggplot object in multiple high-resolution formats (PDF, TIFF, PNG).
#'
#' @param plot The ggplot object to save.
#' @param filename_base The base name for the output files.
#' @param output_dir The directory where plots will be saved.
#' @param width The width of the plot in inches.
#' @param height The height of the plot in inches.
save_publication_plots <- function(plot, filename_base, output_dir, width = 8, height = 6) {
  ggsave(file.path(output_dir, paste0(filename_base, ".pdf")), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename_base, ".tiff")), plot, width = width, height = height, dpi = 300, compression = "lzw")
  ggsave(file.path(output_dir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = 600)
}

# --- Data Loading and Processing ---

network_results_file <- file.path(base_dir, "network_results.rds")
if (!file.exists(network_results_file)) {
  stop("Network results file not found at: ", network_results_file)
}
network_results <- readRDS(network_results_file)

# NOTE: The same network data ('leaf_network' and 'root_network') is used for
# both G1 and G2 genotypes. This will result in identical module structures for
# both genotypes within the same tissue. If genotype-specific networks exist,
# they should be loaded and passed to this function accordingly.
# For example: network_results$leaf_network_g1, network_results$leaf_network_g2
leaf_g1 <- analyze_modules(network_results$leaf_network, "Leaf", "G1")
leaf_g2 <- analyze_modules(network_results$leaf_network, "Leaf", "G2")
root_g1 <- analyze_modules(network_results$root_network, "Root", "G1")
root_g2 <- analyze_modules(network_results$root_network, "Root", "G2")

all_assignments <- bind_rows(leaf_g1, leaf_g2, root_g1, root_g2)

# --- Visualization ---

if (!is.null(all_assignments) && nrow(all_assignments) > 0) {
  module_colors <- colorRampPalette(
    c("#4F6D7A", "#4A86B4", "#5D9C59", "#78BE8A", "#8FC741", "#A9C66E", "#BCCB56", "#D4D156", "#E6E04B", "#F0F06F")
  )(length(unique(all_assignments$module_id)))
  
  module_plot <- ggplot(all_assignments, aes(x = module_size, y = tissue)) +
    geom_point(aes(color = factor(module_id), size = module_size), alpha = 0.8) +
    facet_wrap(~genotype) +
    scale_color_manual(values = module_colors) +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() +
    theme(
      text = element_text(size = 18),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1.5, "lines"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      strip.text = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = "Module-Level Organization in Root & Leaf Networks",
      x = "Module Size",
      y = "Tissue Type",
      color = "Module ID",
      size = "Module Size"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  save_publication_plots(module_plot, "module_organization", output_dir)
  
  cat("Module organization plot saved to:", output_dir, "\n")
  
} else {
  cat("No module assignments to plot.\n")
}

# --- Summary Statistics ---

if (!is.null(all_assignments) && nrow(all_assignments) > 0) {
  summary_file <- file.path(output_dir, "module_summary_statistics.txt")
  
  # Capture summary output and write to a file
  sink(summary_file)
  cat("Summary Statistics for Module Assignments\n")
  cat("=========================================\n\n")
  print(summary(all_assignments))
  sink()
  
  cat("Summary statistics saved to:", summary_file, "\n")
}

cat("Script execution finished.\n")
