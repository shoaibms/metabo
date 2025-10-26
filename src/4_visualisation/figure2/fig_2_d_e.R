# =============================================================================
# Refactored Script: Two-Panel Figure - Panels D & E
# =============================================================================
# Description: Generates a focused two-panel figure:
#   Panel D: Effect-size distributions (Cliff's δ) by tissue and genotype
#   Panel E: Leaf:Root ratio summary (quantified asymmetry)
#
# Author: Plant Systems Biology Lab
# Date: 2025
# =============================================================================

# 1. ENVIRONMENT SETUP
# =============================================================================
# Clear environment and load libraries
rm(list = ls())
gc()

# Load required libraries
required_packages <- c(
  "tidyverse", "ggplot2", "ggridges", "cowplot", "viridis", 
  "RColorBrewer", "scales", "grid", "gridExtra", "broom",
  "effsize", "coin", "boot", "reshape2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 2. FILE PATHS AND DIRECTORIES
# =============================================================================
# Data file paths (matching original color scheme requirements)
data_paths <- list(
  leaf_data = "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_l.csv",
  root_data = "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_r.csv",
  vip_data = "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv",
  temporal_analysis = "C:/Users/ms/Desktop/data_chem_3_10/output/results/initial_stat/initial_stat6e/temporal_analysis.csv"
)

# Output directory
output_dir <- "C:/Users/ms/Desktop/r/chem_data/plot/fig"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Verify all data files exist
missing_files <- data_paths[!sapply(data_paths, file.exists)]
if (length(missing_files) > 0) {
  stop("Missing data files: ", paste(names(missing_files), collapse = ", "))
}

cat("✓ All data files found\n")
cat("✓ Output directory created:", output_dir, "\n")

# 3. DATA LOADING AND PREPROCESSING
# =============================================================================
cat("Loading and preprocessing data...\n")

# Load datasets with robust error handling
safe_read_csv <- function(path, description) {
  tryCatch({
    data <- read.csv(path, stringsAsFactors = FALSE)
    cat("✓ Loaded", description, ":", nrow(data), "rows ×", ncol(data), "columns\n")
    return(data)
  }, error = function(e) {
    stop("Failed to load ", description, ": ", e$message)
  })
}

leaf_data <- safe_read_csv(data_paths$leaf_data, "leaf metabolomics data")
root_data <- safe_read_csv(data_paths$root_data, "root metabolomics data")
vip_data <- safe_read_csv(data_paths$vip_data, "VIP metabolites data")

# Extract metabolite feature columns
get_metabolite_columns <- function(data) {
  grep("^(N_Cluster_|P_Cluster_)", names(data), value = TRUE)
}

leaf_metabolites <- get_metabolite_columns(leaf_data)
root_metabolites <- get_metabolite_columns(root_data)
shared_metabolites <- intersect(leaf_metabolites, root_metabolites)

cat("✓ Found", length(shared_metabolites), "shared metabolites between tissues\n")

# Validate data structure
required_cols <- c("Day", "Genotype", "Treatment")
for (col in required_cols) {
  if (!col %in% names(leaf_data) || !col %in% names(root_data)) {
    stop("Missing required column: ", col)
  }
}

cat("✓ Data validation complete\n")

# 4. PLOTTING CONFIGURATION
# =============================================================================
# Centralized parameters for easy customization of plot aesthetics
plot_config <- list(
  # Color Scheme
  colors = list(
    genotype = c("G1" = "#2ecc71", "G2" = "#33d6d3"),  # Original green and teal
    tissue = c("Leaf" = "#27ae60", "Root" = "#33d6d3"), # Darker green for leaf, teal for root
    accent = "#dee063",  # Red accent for significance
    neutral = c("#34495e", "#7f8c8d", "#bdc3c7")  # Gray palette
  ),
  
  # Font Sizes
  fonts = list(
    base = 12,
    title = 14,
    subtitle = 12,
    tag = 16,
    axis_title = 12,
    axis_text = 11,
    legend_title = 12,
    legend_text = 11,
    strip_text = 12,
    geom_text = 3.5
  ),
  
  # Plot Dimensions (width, height in inches)
  dims = list(
    panel_d = c(10, 4.5),
    panel_e = c(8, 4.5),
    combined = c(12, 4)  # Two panels in one row
  )
)

# 5. VISUAL THEME
# =============================================================================
# Publication theme using centralized parameters
publication_theme <- theme_minimal(base_size = plot_config$fonts$base) +
  theme(
    # Text elements
    plot.title = element_text(size = plot_config$fonts$title, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = plot_config$fonts$subtitle, color = "grey30"),
    plot.tag = element_text(size = plot_config$fonts$tag, face = "bold"),
    plot.tag.position = "topleft",
    axis.title = element_text(size = plot_config$fonts$axis_title, face = "plain"),
    axis.text = element_text(size = plot_config$fonts$axis_text, color = "black"),
    
    # Legend styling
    legend.title = element_text(size = plot_config$fonts$legend_title, face = "bold"),
    legend.text = element_text(size = plot_config$fonts$legend_text),
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),
    
    # Panel and grid
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    
    # Strips for facets
    strip.text = element_text(size = plot_config$fonts$strip_text, face = "plain", color = "black"),
    strip.background = element_rect(fill = "grey95", color = "black"),
    
    # Margins
    plot.margin = margin(15, 15, 15, 15)
  )

# 6. ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate Cliff's δ Effect Sizes
calculate_effect_sizes <- function(leaf_data, root_data, shared_metabolites) {
  cat("Calculating Cliff's δ effect sizes...\n")
  
  results <- list()
  
  for (genotype in c("G1", "G2")) {
    for (tissue_data in list(list(data = leaf_data, name = "Leaf"), 
                             list(data = root_data, name = "Root"))) {
      
      tissue_name <- tissue_data$name
      data <- tissue_data$data
      
      # Get control and treatment data
      control_data <- data %>%
        filter(Genotype == genotype, Treatment == 0) %>%
        select(all_of(shared_metabolites))
      
      treatment_data <- data %>%
        filter(Genotype == genotype, Treatment == 1) %>%
        select(all_of(shared_metabolites))
      
      if (nrow(control_data) > 0 && nrow(treatment_data) > 0) {
        
        for (metabolite in shared_metabolites) {
          control_vals <- control_data[[metabolite]]
          treatment_vals <- treatment_data[[metabolite]]
          
          # Remove NAs
          control_clean <- control_vals[!is.na(control_vals)]
          treatment_clean <- treatment_vals[!is.na(treatment_vals)]
          
          if (length(control_clean) >= 3 && length(treatment_clean) >= 3) {
            # Calculate Cliff's δ
            tryCatch({
              cliff_result <- cliff.delta(treatment_clean, control_clean)
              
              results[[paste(genotype, tissue_name, metabolite, sep = "_")]] <- data.frame(
                Genotype = genotype,
                Tissue = tissue_name,
                Metabolite = metabolite,
                Effect_Size = cliff_result$estimate,
                Magnitude = cliff_result$magnitude,
                CI_Lower = cliff_result$conf.int[1],
                CI_Upper = cliff_result$conf.int[2],
                N_Control = length(control_clean),
                N_Treatment = length(treatment_clean)
              )
            }, error = function(e) {
              # Skip metabolites that cause errors
            })
          }
        }
      }
    }
  }
  
  result_df <- do.call(rbind, results)
  cat("✓ Calculated effect sizes for", nrow(result_df), "metabolite-tissue-genotype combinations\n")
  return(result_df)
}

#' Calculate Leaf:Root Response Ratios
calculate_leafroot_ratios <- function(effect_sizes_data) {
  cat("Calculating Leaf:Root response ratios...\n")
  
  # Define effect size categories
  effect_categories <- effect_sizes_data %>%
    mutate(
      Abs_Effect = abs(Effect_Size),
      Category = case_when(
        Abs_Effect >= 0.474 ~ ">= 0.474",
        Abs_Effect >= 0.33 ~ "0.33-0.474",
        Abs_Effect >= 0.147 ~ "0.147-0.33",
        TRUE ~ "< 0.147"
      )
    ) %>%
    filter(Category %in% c(">= 0.474", "0.33-0.474"))  # Focus on medium and large effects
  
  # Count metabolites by category, tissue, and genotype
  ratio_data <- effect_categories %>%
    group_by(Genotype, Category, Tissue) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Tissue, values_from = Count, values_fill = 0) %>%
    mutate(
      Leaf_Root_Ratio = ifelse(Root == 0, Leaf, Leaf / Root),
      Total_Features = Leaf + Root
    )
  
  cat("✓ Calculated ratios for", nrow(ratio_data), "category-genotype combinations\n")
  return(ratio_data)
}

# 7. EXECUTE ANALYSES
# =============================================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("EXECUTING EFFECT SIZE ANALYSES\n")
cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")

# Run analyses
effect_sizes_data <- calculate_effect_sizes(leaf_data, root_data, shared_metabolites)
leafroot_ratios_data <- calculate_leafroot_ratios(effect_sizes_data)

# 8. PANEL GENERATION
# =============================================================================

#' Panel D: Effect Size Distributions (Ridge Plot)
create_panel_d <- function(effect_sizes_data) {
  cat("Creating Panel D: Effect size distributions...\n")
  
  # Add tissue type abbreviation for plot
  plot_data <- effect_sizes_data %>%
    mutate(Tissue_Type = case_when(
      Tissue == "Leaf" ~ "L",
      Tissue == "Root" ~ "R"
    ))
  
  p <- ggplot(plot_data, aes(x = Effect_Size, y = Tissue_Type, fill = Genotype, color = Genotype)) +
    geom_density_ridges(alpha = 0.7, scale = 0.8, rel_min_height = 0.01) +
    
    # Add reference lines for effect size thresholds
    geom_vline(xintercept = c(-0.474, -0.33, -0.147, 0, 0.147, 0.33, 0.474),
               linetype = "dashed", color = "gray60", alpha = 0.7) +
    
    # Styling with original colors
    scale_fill_manual(values = plot_config$colors$genotype, name = "Genotype", aesthetics = c("fill", "color")) +
    scale_y_discrete(labels = c("L" = "Leaf", "R" = "Root"), expand = expansion(mult = c(0.05, 0.35))) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    coord_cartesian(ylim = c(0.7, 2.5)) +
    
    labs(
      tag = "D",
      title = "Effect Size Distribution by Tissue and Genotype",
      subtitle = "Cliff's δ distributions showing tissue-specific response patterns",
      x = "Effect Size (Cliff's δ)",
      y = "Tissue Type"
    ) +
    
    publication_theme +
    theme(legend.position = "top")
  
  return(p)
}

#' Panel E: Leaf:Root Ratio Summary
create_panel_e <- function(leafroot_ratios_data) {
  cat("Creating Panel E: Leaf:Root ratio summary...\n")
  
  p <- ggplot(leafroot_ratios_data, aes(x = Category, y = Leaf_Root_Ratio, fill = Genotype)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    
    # Add reference line at ratio = 1 (equal leaf:root response)
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 1) +
    
    # Add value labels on bars
    geom_text(aes(label = round(Leaf_Root_Ratio, 1)), 
              position = position_dodge(0.8), vjust = -0.3, size = plot_config$fonts$geom_text, fontface = "plain") +
    
    # Styling with original colors
    scale_fill_manual(values = plot_config$colors$genotype, name = "Genotype") +
    scale_y_continuous(limits = c(0, max(leafroot_ratios_data$Leaf_Root_Ratio) * 1.1)) +
    
    labs(
      tag = "E",
      title = "Leaf:Root Response Ratio",
      subtitle = "Quantified tissue asymmetry in metabolic responses",
      x = "Effect Size Category (|δ|)",
      y = "Leaf:Root Ratio"
    ) +
    
    publication_theme +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# 9. GENERATE PANELS
# =============================================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("GENERATING TWO-PANEL FIGURE\n")
cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")

panel_d <- create_panel_d(effect_sizes_data)
panel_e <- create_panel_e(leafroot_ratios_data)

# 10. SAVE INDIVIDUAL PANELS (PNG FORMAT)
# =============================================================================
cat("Saving individual panels as PNG...\n")

# Save function with PNG output
safe_ggsave_png <- function(plot, filename, width, height) {
  tryCatch({
    ggsave(
      filename = file.path(output_dir, paste0(filename, ".png")),
      plot = plot,
      width = width, height = height,
      dpi = 300, units = "in",
      bg = "white"
    )
    cat("✓ Saved:", filename, ".png\n")
  }, error = function(e) {
    cat("✗ Failed to save", filename, ":", e$message, "\n")
  })
}

# Save individual panels
safe_ggsave_png(panel_d, "Fig2_Panel_D_EffectSizeDistribution", plot_config$dims$panel_d[1], plot_config$dims$panel_d[2])
safe_ggsave_png(panel_e, "Fig2_Panel_E_LeafRootRatio", plot_config$dims$panel_e[1], plot_config$dims$panel_e[2])

# 11. CREATE COMBINED TWO-PANEL FIGURE
# =============================================================================
cat("Creating combined two-panel figure...\n")

# Combine panels in one row
combined_figure <- plot_grid(
  panel_d, panel_e,
  ncol = 2, nrow = 1,
  align = 'h',
  rel_widths = c(1.2, 1)  # Panel D slightly wider for ridge plot
)

# Save combined figure
safe_ggsave_png(combined_figure, "Fig2_TwoPanel_DE_Combined", plot_config$dims$combined[1], plot_config$dims$combined[2])

# 12. SAVE PROCESSED DATA
# =============================================================================
cat("Saving processed data...\n")

write.csv(effect_sizes_data, file.path(output_dir, "processed_effect_sizes_twopanel.csv"), row.names = FALSE)
write.csv(leafroot_ratios_data, file.path(output_dir, "processed_leafroot_ratios_twopanel.csv"), row.names = FALSE)

# 13. SUMMARY STATISTICS
# =============================================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("TWO-PANEL FIGURE ANALYSIS SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")

# Key findings summary
large_effects <- effect_sizes_data %>% filter(abs(Effect_Size) >= 0.474)
medium_effects <- effect_sizes_data %>% filter(abs(Effect_Size) >= 0.33 & abs(Effect_Size) < 0.474)

g1_large_leaf <- sum(large_effects$Genotype == "G1" & large_effects$Tissue == "Leaf")
g1_large_root <- sum(large_effects$Genotype == "G1" & large_effects$Tissue == "Root")

cat("Key Findings:\n")
cat("• Total effect sizes calculated:", nrow(effect_sizes_data), "\n")
cat("• Large effects (|δ| >= 0.474):", nrow(large_effects), "\n")
cat("• Medium effects (0.33 <= |δ| < 0.474):", nrow(medium_effects), "\n")
cat("• G1 large effects - Leaf:", g1_large_leaf, "Root:", g1_large_root, "\n")
if (g1_large_root > 0) {
  cat("• G1 Leaf:Root ratio (large effects):", round(g1_large_leaf / g1_large_root, 1), "\n")
}
cat("• Ratio categories analyzed:", nrow(leafroot_ratios_data), "\n")
cat("• Output files generated: 5\n")

cat("\n✅ Two-Panel Figure Script completed successfully!\n")
cat("📁 All outputs saved to:", output_dir, "\n")
cat("🎯 Combined figure shows Panel D and Panel E in one horizontal row\n")