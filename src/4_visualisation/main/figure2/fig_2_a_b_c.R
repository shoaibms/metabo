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
#
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

#' Calculate Cross-tissue Correlation with Bootstrap CI (Fixed Metabolite Set)
calculate_correlation_with_ci <- function(resilience_data, n_boot = 5000, min_fixed_metabolites = 10) {
  cat("Calculating cross-tissue correlations with fixed metabolite set and Bootstrap CIs...\n")
  set.seed(20250816)
  
  required_tissues <- c("L", "R")
  available_tissues <- sort(unique(as.character(resilience_data$Tissue)))
  if (!all(required_tissues %in% available_tissues)) {
    stop("Expected Tissue labels 'L' and 'R' were not both found in resilience_index.csv")
  }
  
  genotypes <- sort(unique(as.character(resilience_data$Genotype)))
  time_points <- sort(unique(resilience_data$Day))
  
  # A fixed metabolite set is required so temporal rho values are directly comparable.
  # Preferred: one metabolite vector shared across both genotypes, both tissues, and all time points.
  expected_global_cells <- length(genotypes) * length(time_points) * length(required_tissues)
  global_fixed_metabolites <- resilience_data %>%
    filter(Tissue %in% required_tissues, !is.na(Resilience_Index)) %>%
    distinct(Genotype, Day, Tissue, Metabolite) %>%
    group_by(Metabolite) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    filter(n_cells == expected_global_cells) %>%
    arrange(Metabolite) %>%
    pull(Metabolite)
  
  use_global_fixed_set <- length(global_fixed_metabolites) >= min_fixed_metabolites
  if (use_global_fixed_set) {
    cat("✓ Using one fixed metabolite set across both genotypes:", length(global_fixed_metabolites), "metabolites\n")
  } else {
    cat("⚠ Global fixed set has", length(global_fixed_metabolites), "metabolites; using genotype-specific fixed sets instead\n")
  }
  
  # Define the function that boot will use to calculate correlation on a resampled metabolite set.
  # Some bootstrap resamples can be uninformative if one vector has no rank variation;
  # return NA for those draws and remove them when computing percentile CIs.
  boot_correlation_stat <- function(data, indices) {
    sample_data <- data[indices, ]
    if (length(unique(sample_data$L)) < 2 || length(unique(sample_data$R)) < 2) return(NA_real_)
    suppressWarnings(cor(sample_data$L, sample_data$R, method = "spearman", use = "complete.obs"))
  }
  
  results_list <- list()
  diagnostics_list <- list()
  
  for (g in genotypes) {
    if (use_global_fixed_set) {
      fixed_metabolites <- global_fixed_metabolites
      fixed_scope <- "global_across_genotypes"
    } else {
      expected_genotype_cells <- length(time_points) * length(required_tissues)
      fixed_metabolites <- resilience_data %>%
        filter(Genotype == g, Tissue %in% required_tissues, !is.na(Resilience_Index)) %>%
        distinct(Day, Tissue, Metabolite) %>%
        group_by(Metabolite) %>%
        summarise(n_cells = n(), .groups = "drop") %>%
        filter(n_cells == expected_genotype_cells) %>%
        arrange(Metabolite) %>%
        pull(Metabolite)
      fixed_scope <- "genotype_specific"
      cat("  ", g, ": using genotype-specific fixed set of", length(fixed_metabolites), "metabolites\n")
    }
    
    if (length(fixed_metabolites) < min_fixed_metabolites) {
      warning(paste0("Fixed metabolite set for ", g, " has only ", length(fixed_metabolites),
                     " metabolites; correlation and bootstrap CI may be unstable."))
    }
    
    for (tp in time_points) {
      group_data <- resilience_data %>%
        filter(Genotype == g, Day == tp, Metabolite %in% fixed_metabolites, Tissue %in% required_tissues) %>%
        select(Metabolite, Tissue, Resilience_Index) %>%
        pivot_wider(names_from = Tissue, values_from = Resilience_Index) %>%
        arrange(match(Metabolite, fixed_metabolites)) %>%
        drop_na(L, R)
      
      diagnostics_list[[paste(g, tp, sep = "_")]] <- data.frame(
        Genotype = g,
        Time_Point = tp,
        Fixed_Set_Scope = fixed_scope,
        N_Fixed_Metabolites = length(fixed_metabolites),
        N_Used_Metabolites = nrow(group_data)
      )
      
      if (nrow(group_data) > min_fixed_metabolites) {
        rho_observed <- suppressWarnings(cor(group_data$L, group_data$R, method = "spearman", use = "complete.obs"))
        
        boot_results <- boot(
          data = group_data,
          statistic = boot_correlation_stat,
          R = n_boot
        )
        
        boot_values <- as.numeric(boot_results$t[, 1])
        boot_values <- boot_values[is.finite(boot_values)]
        n_valid_boot <- length(boot_values)
        
        if (n_valid_boot >= 100) {
          ci_values <- as.numeric(quantile(boot_values, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
        } else {
          ci_values <- c(NA_real_, NA_real_)
          warning(paste0("Only ", n_valid_boot, " valid bootstrap values for ", g,
                         " time point ", tp, "; CI set to NA."))
        }
        
        results_list[[paste(g, tp, sep = "_")]] <- data.frame(
          Genotype = g,
          Day = tp,
          Correlation = rho_observed,
          CI_Lower = ci_values[1],
          CI_Upper = ci_values[2],
          N_Common_Metabolites = nrow(group_data),
          N_Valid_Bootstrap = n_valid_boot,
          Fixed_Set_Scope = fixed_scope
        )
      } else {
        warning(paste0("Skipping ", g, " time point ", tp,
                       ": only ", nrow(group_data), " matched metabolites after fixed-set filtering."))
      }
    }
  }
  
  final_results <- do.call(rbind, results_list)
  fixed_set_diagnostics <- do.call(rbind, diagnostics_list)
  attr(final_results, "fixed_set_diagnostics") <- fixed_set_diagnostics
  
  cat("✓ Fixed-set correlations with Bootstrap CIs calculated.\n")
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
  
  # Prepare data for the specific G1 callout using first and last sampled time point
  g1_series <- correlation_results_with_ci %>%
    filter(Genotype == "G1") %>%
    arrange(Day)
  g1_vals <- data.frame(
    start_tp = min(g1_series$Day),
    end_tp = max(g1_series$Day),
    start_rho = g1_series$Correlation[which.min(g1_series$Day)],
    end_rho = g1_series$Correlation[which.max(g1_series$Day)]
  )
  
  # Construct the precise label string
  callout_label <- paste0("G1: ρ ", round(g1_vals$start_rho, 3), "→", round(g1_vals$end_rho, 3), " (95% CI)")
  
  panel_a_ci_diagnostics <- correlation_results_with_ci %>%
    mutate(
      CI_Finite = is.finite(CI_Lower) & is.finite(CI_Upper),
      CI_Lower_Below_Display = CI_Lower < 0.20,
      CI_Upper_Above_Display = CI_Upper > 0.60,
      CI_Width = CI_Upper - CI_Lower
    )
  
  write.csv(panel_a_ci_diagnostics,
            file.path(output_dir, "panel_a_ci_diagnostics.csv"),
            row.names = FALSE)
  
  cat("Panel A CI diagnostics written to:", file.path(output_dir, "panel_a_ci_diagnostics.csv"), "\n")
  print(panel_a_ci_diagnostics %>% select(Genotype, Day, Correlation, CI_Lower, CI_Upper, CI_Finite,
                                          CI_Lower_Below_Display, CI_Upper_Above_Display,
                                          N_Common_Metabolites, N_Valid_Bootstrap))
  
  p <- ggplot(correlation_results_with_ci, aes(x = Day, y = Correlation, color = Genotype, fill = Genotype, group = Genotype)) +
    
    # 95% bootstrap CI ribbon. coord_cartesian() below clips visually without dropping
    # rows whose CI extends outside the displayed y-range; using scale limits would
    # silently remove those rows and make ribbons disappear.
    geom_ribbon(
      data = correlation_results_with_ci %>% filter(is.finite(CI_Lower), is.finite(CI_Upper)),
      aes(ymin = CI_Lower, ymax = CI_Upper),
      alpha = 0.12,
      color = NA
    ) +
    
    # Main correlation lines and points with specified line width
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    
    # Add the precise callout for G1 decoupling trend using annotate
    # Arrow coordinates are translated from the matplotlib example
    annotate(
      "segment",
      x = g1_vals$start_tp + 0.1, xend = g1_vals$end_tp,
      y = 0.57, yend = g1_vals$end_rho + 0.015,
      arrow = arrow(length = unit(0.2, "cm"), type = "open", ends="last"),
      color = "black",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = g1_vals$start_tp + 0.05, y = 0.575, # Adjusted position slightly for better alignment
      label = callout_label,
      hjust = 0,
      vjust = 0,
      color = "black",
      size = 3.5
    ) +
    
    # Apply the specified color palette. Keep the legend compact so it does not
    # overlap with Panel B in the horizontal layout.
    scale_color_manual(values = nature_colors, name = "Genotype") +
    scale_fill_manual(values = nature_colors, guide = "none") +
    
    # Set breaks for clarity. Use coord_cartesian rather than scale limits so
    # CI ribbons are not dropped when intervals extend beyond the displayed range.
    scale_x_continuous(breaks = 1:3) +
    scale_y_continuous(
      breaks = seq(0.2, 0.6, 0.1),
      name = "Cross-tissue Correlation (Spearman ρ)"
    ) +
    coord_cartesian(ylim = c(0.15, 0.62)) +
    
    # Final titles and theming
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Time point"
    ) +
    nature_theme +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 8, face = "plain"),
      legend.text = element_text(size = 8),
      legend.key.width = unit(0.45, "cm"),
      legend.spacing.x = unit(0.08, "cm"),
      legend.margin = margin(t = 3),
      plot.margin = margin(5, 10, 5, 5)
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
      title = NULL,
      subtitle = NULL,
      x = "Time point",
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
      title = NULL,
      subtitle = NULL,
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
  labels = c("A", "B", "C"),
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


# Save the dynamically calculated fixed-set correlation data with CIs
write.csv(correlation_results_with_ci, file.path(output_dir, "thesis_correlation_data_with_CIs.csv"), row.names = FALSE)

# Save fixed-set diagnostics for transparent response to Reviewer 2
fixed_set_diagnostics <- attr(correlation_results_with_ci, "fixed_set_diagnostics")
if (!is.null(fixed_set_diagnostics)) {
  write.csv(fixed_set_diagnostics, file.path(output_dir, "fixed_metabolite_set_diagnostics.csv"), row.names = FALSE)
}

cat("\n✅ High-impact Figure 2 script completed successfully!\n")
cat("📁 Outputs saved to:", output_dir, "\n")
cat("📊 Layout: 3 panels in a row, Panel A wider, tissues stacked vertically in B & C\n")
