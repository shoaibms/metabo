# ---
# title: "Multi-level Network Validation and Visualization"
# author: "Your Name/Lab Name"
# date: "Last updated: 2023-10-27"
#
# description: >
#   This script performs a multi-level validation of metabolic networks built from
#   temporal data. It assesses three key aspects of network quality:
#   1. Module Preservation: Evaluates the conservation of network modules across
#      different conditions.
#   2. Network Evolution: Tracks the stability and hub persistence of networks
#      over time.
#   3. Statistical Cross-Validation: Compares effect sizes with permutation-based
#      scores to ensure statistical robustness.
#
#   The script generates a composite figure (Figure 3) summarizing these validation
#   metrics and exports the underlying data for reproducibility.
#
# input:
#   - A single .rds file containing the results of a temporal network analysis.
#     This object should include edge probabilities, null distributions, and
#     permutation results for different tissues and genotypes.
#
# output:
#   - "Fig3_network_validation.pdf": High-resolution composite plot for publication.
#   - "Fig3_network_validation.png": High-resolution composite plot for presentations.
#   - "module_level_validation.csv": Data for panel A.
#   - "network_evolution.csv": Data for panel B.
#   - "statistical_cross_validation.csv": Data for panel C.
#   - "validation_summary_metrics.csv": High-level summary of validation scores.
# ---

# 1. SETUP
# -----------------------------------------------------------------------------
# Load required libraries
library(tidyverse)
library(igraph)
library(boot)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(viridis)

# --- Configuration ---

# NOTE: Hard-coded paths are retained as per user request.
# For greater portability, consider using the 'here' package.
BASE_DIR <- "C:/Users/ms/Desktop/r/chem_data/final"
NETWORK_RESULTS_PATH <- file.path(BASE_DIR, "baysian_new_crosstalk5_V5_5000/temporal_analysis.rds")
OUTPUT_DIR <- file.path(BASE_DIR, "validation_results")

# Create the output directory if it doesn't exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Plotting Theme and Colors ---
NATURE_THEME <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 12, face = "bold")
  )

TISSUE_COLORS <- c("Leaf" = "#2ecc71", "Root" = "#33d6d3")
GENOTYPE_COLORS <- c(
  "Leaf.G1" = "#2ecc71", "Leaf.G2" = "#27ae60",
  "Root.G1" = "#33d6d3", "Root.G2" = "#1e597d"
)

# --- Analysis Constants ---
GENOTYPE_THRESHOLDS <- c("G1" = 0.5, "G2" = 0.4)
PRESERVATION_WEIGHTS <- c(density = 0.4, transitivity = 0.3, modularity = 0.3)
N_GENOTYPE_SAMPLES <- 250
N_VALIDATION_SAMPLES <- 100
TIME_POINTS <- 1:3


# 2. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

#' Calculate Module Preservation Score
#'
#' @param edge_prob A matrix of edge probabilities.
#' @param genotype A string indicating the genotype ("G1" or "G2").
#' @return A single numeric preservation score.
calculate_preservation <- function(edge_prob, genotype) {
  threshold <- GENOTYPE_THRESHOLDS[genotype]
  adj_matrix <- edge_prob > threshold
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

  density <- edge_density(graph)
  transitivity <- transitivity(graph, type = "global")
  modularity <- modularity(cluster_fast_greedy(graph))

  score <- (PRESERVATION_WEIGHTS["density"] * density +
            PRESERVATION_WEIGHTS["transitivity"] * transitivity +
            PRESERVATION_WEIGHTS["modularity"] * modularity)

  return(score)
}


# 3. CORE ANALYSIS FUNCTIONS
# -----------------------------------------------------------------------------

#' Calculate Multi-level Network Validation Metrics
#'
#' @param network_results An object containing all network analysis results.
#' @return A list of data frames with validation metrics.
calculate_all_validation_metrics <- function(network_results) {
  # --- a) Module-Level Statistics ---
  g1_indices <- 1:N_GENOTYPE_SAMPLES
  g2_indices <- (N_GENOTYPE_SAMPLES + 1):(2 * N_GENOTYPE_SAMPLES)

  module_stats <- expand.grid(
    Tissue = c("Leaf", "Root"),
    Genotype = c("G1", "G2")
  ) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      Preservation = {
        edge_prob <- if (Tissue == "Leaf") network_results$leaf_network$edge_prob else network_results$root_network$edge_prob
        calculate_preservation(edge_prob, Genotype)
      },
      Integration = {
        dist <- if (Tissue == "Leaf") network_results$leaf_network$null_distribution else network_results$root_network$null_distribution
        indices <- if (Genotype == "G1") g1_indices else g2_indices
        mean(dist$effect_sizes[indices])
      },
      Robustness = {
        dist <- if (Tissue == "Leaf") network_results$leaf_network$null_distribution else network_results$root_network$null_distribution
        indices <- if (Genotype == "G1") g1_indices else g2_indices
        1 - mean(dist$p_values[indices])
      }
    ) %>%
    ungroup()

  # --- b) Stability and Persistence Metrics ---
  stability_metrics <- expand.grid(
    Time = TIME_POINTS,
    Tissue = c("Leaf", "Root"),
    Genotype = c("G1", "G2")
  ) %>%
    as_tibble() %>%
    arrange(Tissue, Genotype, Time)

  # Extract and assign stability scores in the correct order
  stability_scores <- c(
    network_results$permutation_results$leaf_null$statistics[1:3], # Leaf G1
    network_results$permutation_results$leaf_null$statistics[4:6], # Leaf G2
    network_results$permutation_results$root_null$statistics[1:3], # Root G1
    network_results$permutation_results$root_null$statistics[4:6]  # Root G2
  )

  hub_persistence_scores <- c(
    1 - network_results$permutation_results$leaf_null$p_values[1:3], # Leaf G1
    1 - network_results$permutation_results$leaf_null$p_values[4:6], # Leaf G2
    1 - network_results$permutation_results$root_null$p_values[1:3], # Root G1
    1 - network_results$permutation_results$root_null$p_values[4:6]  # Root G2
  )

  stability_metrics <- stability_metrics %>%
    mutate(
      Stability = stability_scores,
      Hub_Persistence = hub_persistence_scores
    )

  # --- c) Statistical Validation Metrics ---
  validation_stats <- tibble(
    Tissue = rep(c("Leaf", "Root"), each = N_VALIDATION_SAMPLES),
    Effect_Size = c(
      network_results$leaf_network$null_distribution$effect_sizes[1:N_VALIDATION_SAMPLES],
      network_results$root_network$null_distribution$effect_sizes[1:N_VALIDATION_SAMPLES]
    ),
    Permutation_Score = c(
      1 - network_results$leaf_network$null_distribution$p_values[1:N_VALIDATION_SAMPLES],
      1 - network_results$root_network$null_distribution$p_values[1:N_VALIDATION_SAMPLES]
    )
  )

  return(list(
    module_stats = module_stats,
    stability_metrics = stability_metrics,
    validation_stats = validation_stats
  ))
}

#' Export Validation Data to CSV Files
#'
#' @param validation_metrics A list of data frames from the main calculation function.
#' @param output_dir Path to the directory where files will be saved.
export_validation_data <- function(validation_metrics, output_dir) {
  # Panel A data
  write.csv(
    validation_metrics$module_stats %>%
      select(
        Tissue, Genotype,
        Module_Preservation = Preservation,
        Integration_Score = Integration,
        Robustness
      ),
    file = file.path(output_dir, "module_level_validation.csv"),
    row.names = FALSE
  )

  # Panel B data
  write.csv(
    validation_metrics$stability_metrics %>%
      select(
        Time, Tissue, Genotype,
        Network_Stability = Stability,
        Hub_Persistence
      ),
    file = file.path(output_dir, "network_evolution.csv"),
    row.names = FALSE
  )

  # Panel C data
  write.csv(
    validation_metrics$validation_stats,
    file = file.path(output_dir, "statistical_cross_validation.csv"),
    row.names = FALSE
  )
}


# 4. VISUALIZATION FUNCTIONS
# -----------------------------------------------------------------------------

#' Create Panel A: Module-Level Validation Plot
create_module_validation_plot <- function(module_stats) {
  ggplot(module_stats, aes(x = Preservation, y = Integration, size = Robustness)) +
    geom_point(aes(color = Tissue, shape = Genotype), alpha = 0.8, stroke = 1.2) +
    scale_color_manual(values = TISSUE_COLORS) +
    scale_shape_manual(values = c("G1" = 16, "G2" = 17)) +
    scale_size_continuous(range = c(4, 12)) +
    NATURE_THEME +
    labs(
      title = "a  Module-Level Network Validation",
      x = "Module Preservation Score",
      y = "Integration Score",
      color = "Tissue",
      shape = "Genotype",
      size = "Robustness"
    )
}

#' Create Panel B: Network Evolution Analysis Plot
create_network_evolution_plot <- function(stability_metrics) {
  ggplot(stability_metrics, aes(x = Time, y = Stability, group = interaction(Tissue, Genotype), color = interaction(Tissue, Genotype))) +
    geom_line(linewidth = 1.2) +
    geom_point(aes(size = Hub_Persistence), fill = "white", shape = 21, stroke = 1.2) +
    scale_color_manual(values = GENOTYPE_COLORS, name = "Tissue & Genotype") +
    scale_size_continuous(range = c(3, 8), name = "Hub Persistence") +
    NATURE_THEME +
    labs(
      title = "b  Network Evolution Analysis",
      x = "Time Point",
      y = "Network Stability Score"
    )
}

#' Create Panel C: Statistical Cross-Validation Plot
create_cross_validation_plot <- function(validation_stats) {
  ggplot(validation_stats, aes(x = Effect_Size, y = Permutation_Score, color = Tissue)) +
    geom_point(alpha = 0.5) +
    geom_density_2d(color = "white", alpha = 0.5) +
    scale_color_manual(values = TISSUE_COLORS) +
    NATURE_THEME +
    labs(
      title = "c  Statistical Cross-Validation",
      x = "Effect Size",
      y = "Permutation Score"
    )
}


# 5. MAIN EXECUTION
# -----------------------------------------------------------------------------
main <- function() {
  # Load pre-computed network analysis results
  network_results <- readRDS(NETWORK_RESULTS_PATH)

  # Calculate all validation metrics
  validation_metrics <- calculate_all_validation_metrics(network_results)

  # Export detailed data tables for reproducibility
  export_validation_data(validation_metrics, OUTPUT_DIR)

  # Generate individual plot panels
  p1 <- create_module_validation_plot(validation_metrics$module_stats)
  p2 <- create_network_evolution_plot(validation_metrics$stability_metrics)
  p3 <- create_cross_validation_plot(validation_metrics$validation_stats)

  # Combine plots into a single figure
  combined_plot <- plot_grid(
    p1, p2, p3,
    ncol = 3,
    align = 'h',
    axis = 'tblr'
  )

  # Add a publication-quality title to the composite figure
  figure_title <- ggdraw() +
    draw_label(
      "Figure 3 | Multi-level Network Validation Framework",
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(plot.margin = margin(0, 0, 10, 7))

  final_plot <- plot_grid(
    figure_title, combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )

  # Save the final figure in high-resolution formats
  ggsave(
    file.path(OUTPUT_DIR, "Fig3_network_validation.pdf"),
    final_plot, width = 15, height = 5, device = cairo_pdf
  )
  ggsave(
    file.path(OUTPUT_DIR, "Fig3_network_validation.png"),
    final_plot, width = 15, height = 5, dpi = 300, type = "cairo"
  )

  # Save a summary of key validation metrics
  summary_metrics <- tibble(
    Metric = c("Mean Module Preservation", "Mean Network Stability", "Mean Effect Size"),
    Value = c(
      mean(validation_metrics$module_stats$Preservation),
      mean(validation_metrics$stability_metrics$Stability),
      mean(validation_metrics$validation_stats$Effect_Size)
    )
  )

  write.csv(
    summary_metrics,
    file = file.path(OUTPUT_DIR, "validation_summary_metrics.csv"),
    row.names = FALSE
  )

  cat("Analysis complete. Outputs saved to:", OUTPUT_DIR, "\n")
}

# Execute the main analysis pipeline
main()
