# =============================================================================
# Script 1 (Refactored Version): High-Impact Temporal Dynamics Figure
# =============================================================================
# Description: Generates a cohesive, 3-panel figure for publication using
#              actual analysis data arranged horizontally with vertical
#              tissue stacking in panels B and C.
#
#   Panel A: Cross-tissue correlation of the Resilience Index (with Bootstrap CIs) - WIDER
#   Panel B: Underlying median treatment response trajectories - TISSUES STACKED VERTICALLY
#   Panel C: Distribution of the Resilience Index - TISSUES STACKED VERTICALLY
#
# Author: Plant Systems Biology Lab (Enhanced by AI)
# Date: 2025
# =============================================================================

# 1. ENVIRONMENT SETUP
# =============================================================================
# Clear environment and load libraries
rm(list = ls())
gc()

# Load required libraries
required_packages <- c("tidyverse", "ggplot2", "cowplot", "viridis", "RColorBrewer", "scales", "ggridges", "boot")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 2. FILE PATHS AND DIRECTORIES
# =============================================================================
data_paths <- list(
  resilience_data = "C:/Users/ms/Desktop/data_chem_3_10/output/results/initial_stat/initial_stat6e/resilience_index.csv"
)

# Output directory
output_dir <- "C:/Users/ms/Desktop/r/chem_data/plot/fig/high_impact_figure_final"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Verify data file exists
if (!file.exists(data_paths$resilience_data)) {
  stop("Critical data file not found: resilience_index.csv")
}
cat("✓ All data files found\n")

# 3. DATA LOADING AND PREPROCESSING
# =============================================================================
cat("Loading and preprocessing data...\n")
resilience_data <- read.csv(data_paths$resilience_data)
cat("✓ Loaded resilience data:", nrow(resilience_data), "rows ×", ncol(resilience_data), "columns\n")

# Data validation and cleaning
required_cols <- c("Tissue", "Genotype", "Day", "Metabolite", "Resilience_Index", "Control_median", "Treated_median")
if (!all(required_cols %in% names(resilience_data))) {
  stop("Missing one or more required columns in resilience_index.csv")
}
resilience_data <- resilience_data %>%
  mutate(Genotype = as.factor(Genotype),
         Day = as.numeric(Day),
         Tissue = as.factor(Tissue))
cat("✓ Data validation complete\n")

# 4. VISUAL THEME AND STYLING
# =============================================================================
nature_theme <- theme_minimal(base_size = 10) +  # Reduced base size for horizontal layout
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),  # Reduced from 14
    plot.subtitle = element_text(size = 9, color = "grey30", margin = margin(b = 8)),  # Reduced from 11
    axis.title = element_text(size = 10),  # Reduced from 12
    axis.text = element_text(size = 8, color = "black"),  # Reduced from 10
    legend.title = element_text(size = 9),  # Reduced from 11
    legend.text = element_text(size = 8),  # Reduced from 10
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 10, color = "black"),  # Reduced from 12
    strip.background = element_rect(fill = "grey95", color = "black"),
    plot.margin = margin(5, 5, 5, 5)  # Reduced margins
  )

# --- MODIFICATION START ---
# UPDATED: Use the specified production-ready color palette
nature_colors <- c("G1" = "#3CB371", "G2" = "#2AA9B0")
# --- MODIFICATION END ---
arrow_color <- "black" # Option to change arrow color easily

# 5. ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate Cross-tissue Correlation with Bootstrap CI (Thesis Method)
calculate_correlation_with_ci <- function(resilience_data) {
  cat("Calculating cross-tissue correlations with Bootstrap CIs...\n")
  
  # Define the function that boot will use to calculate correlation on a sample
  boot_correlation_stat <- function(data, indices) {
    sample_data <- data[indices, ]
    return(cor(sample_data$L, sample_data$R, method = "spearman", use = "complete.obs"))
  }
  
  results_list <- list()
  
  # Loop through each genotype and day
  for (g in unique(resilience_data$Genotype)) {
    for (d in unique(resilience_data$Day)) {
      
      # Prepare data for this specific group (common metabolites only)
      group_data <- resilience_data %>%
        filter(Genotype == g, Day == d) %>%
        select(Metabolite, Tissue, Resilience_Index) %>%
        pivot_wider(names_from = Tissue, values_from = Resilience_Index) %>%
        na.omit()
      
      if (nrow(group_data) > 10) { # Ensure enough data for bootstrapping
        # Perform bootstrap
        boot_results <- boot(
          data = group_data,
          statistic = boot_correlation_stat,
          R = 5000 # Number of bootstrap replicates
        )
        
        # Calculate 95% confidence interval
        ci <- boot.ci(boot_results, type = "perc", conf = 0.95)
        
        # Store results
        results_list[[paste(g, d)]] <- data.frame(
          Genotype = g,
          Day = d,
          Correlation = boot_results$t0, # Original correlation
          CI_Lower = ci$percent[4],     # 2.5th percentile
          CI_Upper = ci$percent[5]      # 97.5th percentile
        )
      }
    }
  }
  
  final_results <- do.call(rbind, results_list)
  
  cat("✓ Correlations with dynamic CIs calculated.\n")
  return(final_results)
}

#' Calculate Median Treatment Trajectories
calculate_treatment_trajectories <- function(resilience_data) {
  cat("Calculating median treatment response trajectories...\n")
  
  trajectories <- resilience_data %>%
    group_by(Tissue, Genotype, Day) %>%
    summarise(
      Control = median(Control_median, na.rm = TRUE),
      Treated = median(Treated_median, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    pivot_longer(cols = c(Control, Treated), names_to = "Condition", values_to = "Median_Abundance")
  
  cat("✓ Treatment trajectories calculated\n")
  return(trajectories)
}

# 6. EXECUTE ANALYSES
# =============================================================================
correlation_results_with_ci <- calculate_correlation_with_ci(resilience_data)
treatment_trajectories <- calculate_treatment_trajectories(resilience_data)

# 7. PANEL GENERATION
# =============================================================================

# --- MODIFICATION START ---
# UPDATED: Panel A function replaced entirely with the corrected version
#' Panel A - Cross-tissue Coordination (Statistically Correct & Production-Ready)
create_panel_a <- function(correlation_results_with_ci) {
  cat("Creating Final Panel A with production-ready styling...\n")
  
  # Prepare data for the specific G1 callout
  g1_vals <- correlation_results_with_ci %>%
    filter(Genotype == "G1") %>%
    summarise(
      start_rho = Correlation[Day == 1],
      end_rho = Correlation[Day == 3]
    )
  
  # Construct the precise label string
  callout_label <- paste0("G1: ρ ", round(g1_vals$start_rho, 3), "→", round(g1_vals$end_rho, 3), " (95% CI)")
  
  p <- ggplot(correlation_results_with_ci, aes(x = Day, y = Correlation, color = Genotype, fill = Genotype)) +
    
    # 95% bootstrap CI ribbon with specified alpha
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.18, linewidth = 0) +
    
    # Main correlation lines and points with specified line width
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    
    # Add the precise callout for G1 decoupling trend using annotate
    # Arrow coordinates are translated from the matplotlib example
    annotate(
      "segment",
      x = 1.1, xend = 3,
      y = 0.57, yend = g1_vals$end_rho + 0.015,
      arrow = arrow(length = unit(0.2, "cm"), type = "open", ends="last"),
      color = "black",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = 1.05, y = 0.575, # Adjusted position slightly for better alignment
      label = callout_label,
      hjust = 0,
      vjust = 0,
      color = "black",
      size = 3.5
    ) +
    
    # Apply the specified color palette and legend title
    scale_color_manual(
      values = nature_colors,
      name = "Spearman ρ (mean ± 95% bootstrap CI; 5,000 resamples)"
    ) +
    scale_fill_manual(values = nature_colors, guide = "none") +
    
    # Set breaks for clarity and fix y-axis to specified manuscript style
    scale_x_continuous(breaks = 1:3) +
    scale_y_continuous(limits = c(0.20, 0.60), name = "Cross-tissue Correlation (Spearman ρ)") +
    
    # Final titles and theming
    labs(
      title = "A. Coordinated Resilience Response Decouples in Tolerant Genotype",
      subtitle = "Spearman correlation of metabolic Resilience Index between tissues",
      x = "Time Point (Days)"
    ) +
    nature_theme +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "plain"), # Legend title as plain text
      legend.margin = margin(t = 5)
    )
  
  return(p)
}
# --- MODIFICATION END ---

#' Panel B - Tissue-Specific Treatment Response (VERTICAL TISSUE STACKING)
create_panel_b <- function(treatment_trajectories) {
  cat("Creating Panel B: Tissue-Specific Treatment Response with vertical stacking...\n")
  
  p <- ggplot(treatment_trajectories, aes(x = Day, y = Median_Abundance, color = Condition, linetype = Genotype)) +
    geom_line(size = 1.1) +
    geom_point(size = 2.5, aes(shape = Condition)) +
    # CHANGED: facet_grid instead of facet_wrap for vertical stacking
    facet_grid(Tissue ~ ., scales = "free_y", labeller = labeller(Tissue = c("L" = "Leaf", "R" = "Root"))) +
    scale_color_manual(values = c("Control" = "#1eaeb3", "Treated" = "#b2bf5c"), name = "Condition") +
    scale_linetype_manual(values = c("G1" = "solid", "G2" = "dashed"), name = "Genotype") +
    scale_shape_manual(values = c("Control" = 1, "Treated" = 16), name = "Condition") +
    scale_x_continuous(breaks = 1:3) +
    labs(
      title = "B. Tolerant Genotype Mounts\na Stronger, Coordinated Response",
      subtitle = "Median metabolite abundance\nin control vs. treated plants",
      x = "Time Point (Days)",
      y = "Median Abundance (a.u.)"
    ) +
    nature_theme +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 5)
    )
  
  return(p)
}

#' Panel C - Distribution of the Resilience Index (VERTICAL TISSUE STACKING)
create_panel_c <- function(resilience_data) {
  cat("Creating Panel C: Distribution of the Resilience Index with vertical stacking...\n")
  
  p <- ggplot(resilience_data, aes(x = Resilience_Index, y = Genotype, fill = Genotype)) +
    geom_density_ridges(alpha = 0.8, scale = 1.2, rel_min_height = 0.01, linewidth = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
    # CHANGED: facet_grid instead of facet_wrap for vertical stacking
    facet_grid(Tissue ~ ., labeller = labeller(Tissue = c("L" = "Leaf", "R" = "Root"))) +
    scale_fill_manual(values = nature_colors, guide = "none") +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 10), labels = c("0.1", "0.5", "1", "2", "10")) +
    annotation_logticks(sides = "b") +
    labs(
      title = "C. Resilience is More Suppressed\nin Tolerant Tissues",
      subtitle = "Distribution of Resilience Index\n(Treated / Control) across all time points",
      x = "Resilience Index (Log Scale)",
      y = "Genotype"
    ) +
    nature_theme +
    theme(
      axis.text.y = element_text(),
      legend.position = "none"  # Remove legend since colors are self-explanatory
    )
  
  return(p)
}

# 8. GENERATE & ASSEMBLE FIGURE
# =============================================================================
cat("\nGenerating and assembling final figure...\n")

panel_a <- create_panel_a(correlation_results_with_ci)
panel_b <- create_panel_b(treatment_trajectories)
panel_c <- create_panel_c(resilience_data)

# CHANGED: Combine panels into a single figure with horizontal layout
final_figure <- plot_grid(
  panel_a,
  panel_b,
  panel_c,
  ncol = 3,  # CHANGED: 3 columns instead of 1
  align = 'v',
  labels = "AUTO", # Use auto-labeling for simplicity
  rel_widths = c(1.5, 1, 1),  # ADDED: Panel A gets 1.5x width of others
  label_size = 12
)

# 9. SAVE OUTPUT
# =============================================================================
cat("Saving final figure and processed data...\n")

# --- MODIFICATION START ---
# UPDATED: Save as both vector (PDF) and high-resolution raster (PNG)
ggsave(
  filename = file.path(output_dir, "Figure2_High_Impact_Temporal_Horizontal.pdf"),
  plot = final_figure,
  width = 11, height = 4.5,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "Figure2_High_Impact_Temporal_Horizontal.png"),
  plot = final_figure,
  width = 11, height = 4.5,
  dpi = 300,
  bg = "white"
)
# --- MODIFICATION END ---


# Save the dynamically calculated correlation data with CIs
write.csv(correlation_results_with_ci, file.path(output_dir, "thesis_correlation_data_with_CIs.csv"), row.names = FALSE)

cat("\n✅ High-impact Figure 2 script completed successfully!\n")
cat("📁 Outputs saved to:", output_dir, "\n")
cat("📊 Layout: 3 panels in a row, Panel A wider, tissues stacked vertically in B & C\n")