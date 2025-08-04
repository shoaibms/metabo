# -----------------------------------------------------------------------------
# Script: temporal_corr.R
#
# Description:
# This script performs a temporal analysis of tissue-specific molecular responses
# in two genotypes (G1 and G2) under osmotic stress. It generates a multi-panel
# figure (Figure 2) that includes:
#   - Cross-tissue correlation over time.
#   - Response magnitude distribution in leaf and root tissues.
#   - Feature redistribution between tissues for each genotype.
#
# Author: [Your Name/Lab Name]
# Date: [Date of last modification]
# -----------------------------------------------------------------------------

# 1. SETUP
# -----------------------------------------------------------------------------
# Load required libraries
library(tidyverse)
library(cowplot)
library(boot)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(grid)

# Define file paths
# NOTE: For reproducibility, users should update these paths to their own file locations.
leaf_data_path <- "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_l.csv"
root_data_path <- "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_r.csv"
vip_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/vip_bonferroni/VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv"
metrics_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/network_metrics_summary.csv"
output_dir <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/temporal_analysis"

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
leaf_data <- read.csv(leaf_data_path)
root_data <- read.csv(root_data_path)
vip_data <- read.csv(vip_path)
metrics_data <- read.csv(metrics_path)

# 2. PLOT THEME
# -----------------------------------------------------------------------------
# A consistent theme for all plots, inspired by Nature journals.
nature_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 20, color = "gray30"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 18),
    legend.position = "right"
  )

# 3. ANALYSIS FUNCTIONS
# -----------------------------------------------------------------------------

#' Calculate Temporal Correlation with Bootstrapped Confidence Intervals
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as x.
#' @param n_bootstrap Number of bootstrap replicates.
#' @return A list containing the Spearman correlation and 95% confidence intervals.
calculate_temporal_correlation <- function(x, y, n_bootstrap = 1000) {
  correlation <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")

  boot_cor <- function(data, indices) {
    x_boot <- data$x[indices]
    y_boot <- data$y[indices]
    return(cor(x_boot, y_boot, method = "spearman", use = "pairwise.complete.obs"))
  }

  boot_data <- data.frame(x = x, y = y)
  boot_results <- boot(boot_data, boot_cor, R = n_bootstrap)
  ci <- boot.ci(boot_results, type = "perc")

  return(list(
    correlation = correlation,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

#' Calculate Temporal Metrics Across Tissues
#'
#' This function computes the cross-tissue correlation for each genotype at each time point.
#' @param leaf_data A data frame with leaf metabolomics data.
#' @param root_data A data frame with root metabolomics data.
#' @return A data frame with temporal correlation metrics.
calculate_temporal_metrics <- function(leaf_data, root_data) {
  results <- list()

  get_molecular_features <- function(data) {
    grep("^(N_Cluster_|P_Cluster_)", names(data), value = TRUE)
  }

  leaf_features <- get_molecular_features(leaf_data)
  root_features <- get_molecular_features(root_data)
  common_features <- intersect(leaf_features, root_features)

  for (genotype in c("G1", "G2")) {
    for (day in unique(leaf_data$Day)) {
      leaf_subset <- leaf_data %>%
        filter(Genotype == genotype, Day == day) %>%
        select(all_of(common_features)) %>%
        as.matrix() %>%
        as.vector()

      root_subset <- root_data %>%
        filter(Genotype == genotype, Day == day) %>%
        select(all_of(common_features)) %>%
        as.matrix() %>%
        as.vector()

      if (length(leaf_subset) == length(root_subset) && length(leaf_subset) > 1) {
        cor_results <- calculate_temporal_correlation(leaf_subset, root_subset)

        results[[paste(genotype, day)]] <- tibble(
          Genotype = genotype,
          Day = day,
          Correlation = cor_results$correlation,
          CI_Lower = cor_results$ci_lower,
          CI_Upper = cor_results$ci_upper
        )
      }
    }
  }

  return(bind_rows(results))
}

#' Calculate Response Magnitudes
#'
#' Computes the mean and standard error of the absolute response for each tissue,
#' genotype, and day.
#' @param leaf_data A data frame with leaf metabolomics data.
#' @param root_data A data frame with root metabolomics data.
#' @return A list of two data frames (leaf and root) with response magnitudes.
calculate_response_magnitudes <- function(leaf_data, root_data) {
  process_magnitude <- function(data, tissue_name) {
    molecular_cols <- grep("^(N_Cluster_|P_Cluster_)", names(data), value = TRUE)
    data %>%
      select(Genotype, Day, any_of(molecular_cols)) %>%
      pivot_longer(
        cols = any_of(molecular_cols),
        names_to = "Feature",
        values_to = "Value"
      ) %>%
      group_by(Genotype, Day) %>%
      summarise(
        mean_response = mean(abs(Value), na.rm = TRUE),
        se_response = sd(abs(Value), na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      ) %>%
      mutate(Tissue = tissue_name)
  }

  leaf_magnitude <- process_magnitude(leaf_data, "Leaf")
  root_magnitude <- process_magnitude(root_data, "Root")

  return(list(leaf = leaf_magnitude, root = root_magnitude))
}


# 4. PLOTTING FUNCTIONS
# -----------------------------------------------------------------------------

#' Panel A: Create Temporal Correlation Plot
#' @param temporal_metrics A data frame from calculate_temporal_metrics().
#' @return A ggplot object.
create_temporal_correlation_plot <- function(temporal_metrics) {
  p <- temporal_metrics %>%
    mutate(
      Day = factor(Day),
      Correlation = as.numeric(Correlation),
      CI_Lower = as.numeric(CI_Lower),
      CI_Upper = as.numeric(CI_Upper)
    ) %>%
    ggplot(aes(x = Day, y = Correlation, color = Genotype, group = Genotype)) +
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Genotype), alpha = 0.2, linetype = 0) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
    scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
    scale_y_continuous(limits = c(0.2, 0.4), breaks = seq(0.2, 0.4, by = 0.05)) +
    nature_theme +
    labs(
      title = "Cross-tissue Correlation",
      x = "Time (Days)",
      y = expression(paste("Cross-tissue Correlation (", rho, ")"))
    )
  return(p)
}

#' Panel B: Create Response Magnitude Comparison Plot
#' @param response_magnitudes A list of data frames from calculate_response_magnitudes().
#' @return A ggplot object.
create_magnitude_comparison_plot <- function(response_magnitudes) {
  p <- ggplot() +
    geom_line(data = response_magnitudes$leaf, aes(x = Day, y = mean_response, color = Genotype, linetype = "Leaf"), size = 1) +
    geom_line(data = response_magnitudes$root, aes(x = Day, y = mean_response, color = Genotype, linetype = "Root"), size = 1) +
    geom_point(data = response_magnitudes$leaf, aes(x = Day, y = mean_response, color = Genotype), size = 3) +
    geom_point(data = response_magnitudes$root, aes(x = Day, y = mean_response, color = Genotype), size = 3) +
    geom_errorbar(data = response_magnitudes$leaf, aes(x = Day, ymin = mean_response - se_response, ymax = mean_response + se_response, color = Genotype), width = 0.2) +
    geom_errorbar(data = response_magnitudes$root, aes(x = Day, ymin = mean_response - se_response, ymax = mean_response + se_response, color = Genotype), width = 0.2) +
    scale_color_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3")) +
    scale_linetype_manual(values = c("Leaf" = "solid", "Root" = "dashed")) +
    nature_theme +
    labs(
      title = "Response Magnitude Distribution",
      x = "Time (Days)",
      y = "Response Magnitude",
      linetype = "Tissue",
      color = "Genotype"
    )
  return(p)
}

#' Panel C/D: Create Feature Distribution Plot
#'
#' This function creates a bar plot showing the distribution of resilient features
#' between leaf and root tissues for each genotype.
#' @return A ggplot object.
create_feature_distribution_plot <- function() {
  # This data appears to be pre-calculated.
  feature_dist <- data.frame(
    Group = factor(c("G1", "G1", "G2", "G2")),
    Tissue = factor(c("Leaf", "Root", "Leaf", "Root")),
    Resilient_Features = c(7.80, 13.00, 16.36, 16.95)
  )

  # Prepare annotation data
  df_annot <- feature_dist %>%
    group_by(Group) %>%
    summarise(
      Differential = abs(diff(Resilient_Features)),
      .groups = 'drop'
    ) %>%
    mutate(
      Label = paste0("Delta~", format(round(Differential, 2), nsmall = 2), "~'%'")
    )

  max_value <- max(feature_dist$Resilient_Features)

  p <- ggplot(feature_dist, aes(x = Tissue, y = Resilient_Features, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_segment(
      data = df_annot,
      aes(x = 1, xend = 2, y = max_value * 1.05, yend = max_value * 1.05),
      arrow = arrow(ends = "both", type = "closed", length = unit(0.1, "cm")),
      color = "black",
      size = 0.5,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = df_annot,
      aes(x = 1.5, y = max_value * 1.08, label = Label),
      parse = TRUE,
      color = "black",
      size = 5,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c("G1" = "#2ecc71", "G2" = "#33d6d3"), name = "Genotype") +
    facet_grid(~ Group, scales = "free_x", space = "free") +
    nature_theme +
    theme(
      strip.text = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = "Feature Redistribution",
      x = NULL,
      y = "Resilient Features (%)"
    )

  return(p)
}


# 5. MAIN EXECUTION
# -----------------------------------------------------------------------------

#' Main function to run the analysis and generate the final plot
#' @return A combined ggplot object.
main <- function() {
  # Perform calculations
  temporal_metrics <- calculate_temporal_metrics(leaf_data, root_data)
  response_magnitudes <- calculate_response_magnitudes(leaf_data, root_data)

  # Generate individual plots
  p1 <- create_temporal_correlation_plot(temporal_metrics)
  p2 <- create_magnitude_comparison_plot(response_magnitudes)
  p4 <- create_feature_distribution_plot() # Note: This function uses hardcoded data.
  
  # The call to create_temporal_trends_plot() was removed as it is not defined in this script.
  # If you have this function, you can un-comment the line below and add p3 to the plot_grid call.
  # p3 <- create_temporal_trends_plot()

  # Combine plots into a grid
  # If p3 is available, change ncol to 2 and add p3 to the plot_grid call.
  combined_plot <- plot_grid(
    p1, p2, p4,
    ncol = 2,
    align = 'v',
    axis = 'tblr',
    labels = c("A", "B", "C"), # Adjusted labels
    label_size = 20
  )

  # Create a title for the entire figure
  title <- ggdraw() +
    draw_label(
      "Figure 2 | Temporal dynamics and coordination of tissue-specific molecular responses under osmotic stress",
      fontface = 'bold',
      x = 0.01,
      hjust = 0,
      size = 16
    ) +
    theme(plot.margin = margin(0, 0, 10, 7))

  # Assemble the final plot with title
  final_plot <- plot_grid(
    title, combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )

  # Save the final plot in both PDF and PNG formats
  ggsave(
    file.path(output_dir, "Fig2_temporal_dynamics.pdf"),
    final_plot, width = 14, height = 10, dpi = 300
  )
  ggsave(
    file.path(output_dir, "Fig2_temporal_dynamics.png"),
    final_plot, width = 14, height = 10, dpi = 300
  )

  return(final_plot)
}

# Execute the main function and display the plot
if (interactive()) {
  main()
}
