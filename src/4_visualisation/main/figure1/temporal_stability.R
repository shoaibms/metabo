#' ---
#' title: "Temporal Stability Analysis of Co-expression Networks"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' ---
#'
#' This script analyzes the temporal stability of metabolic networks. It calculates
#' network metrics (edge density, modularity, and hub degree) for different
#' tissues, genotypes, and timepoints. The results are then visualized as a
#' scatter plot to show the integrated stability.

# 1. SETUP ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(igraph)
library(RColorBrewer)

# -- Constants
# File paths (hardcoded as per user request)
BASE_DIR <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"
OUTPUT_DIR <- file.path(BASE_DIR, "temporal_stability")

# Analysis parameters
COR_THRESHOLD <- 0.7
TISSUES <- c("L", "R")
GENOTYPES <- c("G1", "G2")

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}


# 2. DATA LOADING -----------------------------------------------------------

# Read network results from an RDS file
network_results <- readRDS(file.path(BASE_DIR, "processed_data.rds"))
pivot_data <- network_results$pivot_data
DAYS <- sort(unique(pivot_data$Day))


# 3. ANALYSIS ---------------------------------------------------------------

#' Calculate Network Stability Metrics for a Single Timepoint
#'
#' @description
#' This function takes a subset of data for a specific tissue, genotype, and
#' timepoint, and calculates several network stability metrics.
#'
#' @param data_subset A dataframe with metabolite abundance data.
#' @param tissue A string indicating the tissue type (e.g., "L" for leaf).
#' @param genotype A string indicating the genotype (e.g., "G1").
#' @param cor_threshold A numeric value for the correlation threshold to
#'   define edges.
#'
#' @return A dataframe with stability metrics: `edge_consistency`,
#'   `module_preservation`, and `hub_conservation`. Returns `NULL` on error.
#'
calculate_timepoint_stability <- function(data_subset, tissue, genotype, cor_threshold = 0.7) {
  tryCatch({
    # Create a Spearman correlation matrix
    cor_matrix <- cor(data_subset, method = "spearman", use = "pairwise.complete.obs")

    # Create a graph from the absolute correlation matrix values
    graph <- graph_from_adjacency_matrix(
      abs(cor_matrix) > cor_threshold,
      mode = "undirected",
      weighted = TRUE
    )

    # Calculate network metrics
    edge_density_value <- edge_density(graph)

    modularity_value <- tryCatch({
      modularity(cluster_louvain(graph))
    }, error = function(e) NA)

    hub_value <- if (vcount(graph) > 0) {
      mean(sort(degree(graph), decreasing = TRUE)[1:min(10, vcount(graph))])
    } else {
      NA
    }

    data.frame(
      tissue = tissue,
      genotype = genotype,
      edge_consistency = edge_density_value,
      module_preservation = modularity_value,
      hub_conservation = hub_value
    )

  }, error = function(e) {
    warning(paste("Error in calculation for", tissue, genotype, ":", e$message))
    return(NULL)
  })
}

# --- Main Analysis Loop
# Calculate stability metrics across all timepoints, tissues, and genotypes.
stability_results <- list()

for (tissue in TISSUES) {
  for (genotype in GENOTYPES) {
    metabolite_cols <- colnames(pivot_data)[grep(paste0("^", tissue, "_"), colnames(pivot_data))]

    for (day in DAYS) {
      subset_data <- pivot_data %>%
        filter(Day == day, Genotype == genotype) %>%
        select(all_of(metabolite_cols))

      if (ncol(subset_data) > 0 && nrow(subset_data) > 1) {
        metrics <- calculate_timepoint_stability(subset_data, tissue, genotype, COR_THRESHOLD)
        if (!is.null(metrics)) {
          metrics$day <- day
          stability_results[[length(stability_results) + 1]] <- metrics
        }
      }
    }
  }
}

# Combine list of dataframes into a single dataframe
stability_df <- bind_rows(stability_results)


# 4. VISUALIZATION ----------------------------------------------------------

# Create a grouping variable for plotting
stability_df <- stability_df %>%
  mutate(group = paste(tissue, genotype, sep = "_"))

# --- Plotting Theme and Palette
nature_theme <- theme_minimal(base_family = "sans") +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

color_palette <- c(
  "L_G1" = "#34eb81", "L_G2" = "#1e9651",
  "R_G1" = "#33d6d3", "R_G2" = "#1e597d"
)

group_labels <- c(
  "L_G1" = "Leaf G1", "L_G2" = "Leaf G2",
  "R_G1" = "Root G1", "R_G2" = "Root G2"
)

# --- Create Plot
integrated_stability_plot <- ggplot(
  stability_df,
  aes(
    x = edge_consistency,
    y = module_preservation,
    size = hub_conservation,
    color = group
  )
) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    name = "Tissue-Genotype",
    values = color_palette,
    labels = group_labels
  ) +
  scale_size_continuous(name = "Hub Conservation", range = c(3, 10)) +
  labs(
    x = "Network Density",
    y = "Module Preservation"
  ) +
  nature_theme


# 5. SAVE PLOT ----------------------------------------------------------------

output_file_base <- file.path(OUTPUT_DIR, "integrated_stability_plot")

ggsave(
  paste0(output_file_base, ".pdf"),
  integrated_stability_plot, width = 8, height = 6, device = cairo_pdf
)

ggsave(
  paste0(output_file_base, ".tiff"),
  integrated_stability_plot, width = 8, height = 6, dpi = 300
)

ggsave(
  paste0(output_file_base, ".png"),
  integrated_stability_plot, width = 8, height = 6, dpi = 300
)
