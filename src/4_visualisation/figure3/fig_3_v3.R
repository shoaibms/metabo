# =============================================================================
# Refactored Script: Six-Panel Figure - Panels A-F
# =============================================================================
# Description: Generates a comprehensive six-panel figure:
#   Panel A: Hub Persistence Analysis
#   Panel B: VIP Metabolite Tracking
#   Panel C: Feature Change Counts
#   Panel D: Feature Stability
#   Panel E: Mean Response Validation
#   Panel F: Null Distribution Validation
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
# Data file paths
data_paths <- list(
  leaf_data = "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_l.csv",
  root_data = "C:/Users/ms/Desktop/data_chem_3_10/data/data/n_p_r.csv",
  vip_data = "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv",
  temporal_analysis = "C:/Users/ms/Desktop/data_chem_3_10/output/results/initial_stat/initial_stat6e/temporal_analysis.csv"
)

# Hub files paths
hub_paths <- list(
  leaf_g1 = "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/leaf_g1_hub_metabolites.csv",
  leaf_g2 = "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/leaf_g2_hub_metabolites.csv",
  root_g1 = "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/root_g1_hub_metabolites.csv",
  root_g2 = "C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/root_g2_hub_metabolites.csv"
)

# Output directory
output_dir <- "C:/Users/ms/Desktop/r/chem_data/plot/fig"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Verify all data files exist
all_paths <- c(data_paths, hub_paths)
missing_files <- all_paths[!sapply(all_paths, file.exists)]
if (length(missing_files) > 0) {
  cat("Warning: Missing files:\n")
  for (file in names(missing_files)) {
    cat("  ", file, ":", missing_files[[file]], "\n")
  }
}

cat("✓ Output directory created:", output_dir, "\n")

# 3. DATA LOADING AND PREPROCESSING
# =============================================================================
cat("Loading and preprocessing data...\n")

# Load datasets with robust error handling
safe_read_csv <- function(path, description) {
  if (!file.exists(path)) {
    cat("Warning: File not found:", path, "\n")
    return(NULL)
  }
  tryCatch({
    data <- read.csv(path, stringsAsFactors = FALSE)
    cat("✓ Loaded", description, ":", nrow(data), "rows ×", ncol(data), "columns\n")
    return(data)
  }, error = function(e) {
    cat("Warning: Failed to load", description, ":", e$message, "\n")
    return(NULL)
  })
}

# Load main datasets
leaf_data <- safe_read_csv(data_paths$leaf_data, "leaf metabolomics data")
root_data <- safe_read_csv(data_paths$root_data, "root metabolomics data")
vip_data <- safe_read_csv(data_paths$vip_data, "VIP metabolites data")
temporal_data <- safe_read_csv(data_paths$temporal_analysis, "temporal analysis data")

# Load hub datasets
hub_data <- list()
for (condition in names(hub_paths)) {
  hub_data[[condition]] <- safe_read_csv(hub_paths[[condition]], paste("hub data for", condition))
}

# Remove NULL entries
hub_data <- hub_data[!sapply(hub_data, is.null)]

# Extract metabolite feature columns
get_metabolite_columns <- function(data) {
  if (is.null(data)) return(character(0))
  grep("^(N_Cluster_|P_Cluster_)", names(data), value = TRUE)
}

leaf_metabolites <- get_metabolite_columns(leaf_data)
root_metabolites <- get_metabolite_columns(root_data)
shared_metabolites <- intersect(leaf_metabolites, root_metabolites)

cat("✓ Found", length(shared_metabolites), "shared metabolites between tissues\n")

# 4. PLOTTING CONFIGURATION
# =============================================================================
# Original color palette and plot settings
color_palette <- list(
  genotype = c("G1" = "#2ecc71", "G2" = "#33d6d3"),  # Green and teal
  tissue = c("Leaf" = "#27ae60", "Root" = "#33d6d3"), # Forest green, teal
  hub_persistence = c("Shared" = "#e74c3c", "G1_only" = "#2ecc71", "G2_only" = "#33d6d3"),
  vip_tracking = c("Increasing" = "#e74c3c", "Decreasing" = "#3498db", "Stable" = "#95a5a6"),
  feature_change = c("Significant" = "#e67e22", "Non_significant" = "#bdc3c7"),
  accent = "#dee063",
  neutral = c("#34495e", "#7f8c8d", "#bdc3c7")
)

# Plot Sizing
plot_sizing <- list(
  individual_panel_width = 6,
  individual_panel_height = 4.5,
  combined_panel_width = 16,  # Reduced width for tighter layout
  combined_panel_height = 10,  # Reduced height
  
  base_text_size = 11,
  title_text_size = 13,
  subtitle_text_size = 11,
  axis_title_size = 11,
  axis_text_size = 10,
  legend_title_size = 11,
  legend_text_size = 10,
  strip_text_size = 11,
  geom_text_size = 3.5,
  axis_y_text_size = 9,
  placeholder_text_size = 6,
  combined_label_size = 14
)

# 5. VISUAL THEME
# =============================================================================
create_nature_theme <- function(plot_sizing, y_title_margin = 10) {
  theme_minimal(base_size = plot_sizing$base_text_size) +
    theme(
      # Text elements
      plot.title = element_text(size = plot_sizing$title_text_size, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = plot_sizing$subtitle_text_size, color = "grey30"),
      axis.title = element_text(size = plot_sizing$axis_title_size, face = "plain"),
      axis.title.y = element_text(size = plot_sizing$axis_title_size, face = "plain", margin = margin(r = y_title_margin)),
      axis.text = element_text(size = plot_sizing$axis_text_size, color = "black"),
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = "topleft",
      
      # Legend styling
      legend.title = element_text(size = plot_sizing$legend_title_size, face = "bold"),
      legend.text = element_text(size = plot_sizing$legend_text_size),
      legend.position = "top",  # Changed to top for all panels
      legend.box.spacing = unit(0.3, "cm"),  # Reduced spacing
      
      # Panel and grid
      panel.grid.major = element_line(color = "grey92", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      
      # Strips for facets
      strip.text = element_text(size = plot_sizing$strip_text_size, face = "plain", color = "black"),
      strip.background = element_rect(fill = "grey95", color = "black"),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10)  # Reduced margins
    )
}

# 6. ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate Cliff's Delta Effect Sizes
calculate_effect_sizes <- function(leaf_data, root_data, shared_metabolites) {
  cat("Calculating Cliff's Delta effect sizes...\n")
  
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
            # Calculate Cliff's Delta
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

#' Calculate Hub Persistence Analysis
calculate_hub_persistence <- function(hub_data) {
  cat("Calculating hub persistence analysis...\n")
  
  if (length(hub_data) == 0) {
    cat("Warning: No hub data available\n")
    return(NULL)
  }
  
  # Extract hub metabolites for each condition
  hub_sets <- list()
  for (condition in names(hub_data)) {
    if (!is.null(hub_data[[condition]]) && nrow(hub_data[[condition]]) > 0) {
      if ("Metabolite" %in% names(hub_data[[condition]])) {
        hub_sets[[condition]] <- hub_data[[condition]]$Metabolite
      } else {
        # Try first column if Metabolite column doesn't exist
        hub_sets[[condition]] <- hub_data[[condition]][, 1]
      }
    }
  }
  
  # Create a standardized data frame structure
  results <- data.frame(
    Analysis_Type = character(),
    Condition1 = character(),
    Condition2 = character(),
    Group1_hubs = numeric(),
    Group2_hubs = numeric(),
    Shared_hubs = numeric(),
    Group1_unique = numeric(),
    Group2_unique = numeric(),
    Jaccard_similarity = numeric(),
    Comparison_Label = character(),
    stringsAsFactors = FALSE
  )
  
  # Analyze overlap between genotypes within tissues
  for (tissue in c("leaf", "root")) {
    g1_condition <- paste(tissue, "g1", sep = "_")
    g2_condition <- paste(tissue, "g2", sep = "_")
    
    if (g1_condition %in% names(hub_sets) && g2_condition %in% names(hub_sets)) {
      g1_hubs <- hub_sets[[g1_condition]]
      g2_hubs <- hub_sets[[g2_condition]]
      
      if (length(g1_hubs) > 0 && length(g2_hubs) > 0) {
        shared_hubs <- intersect(g1_hubs, g2_hubs)
        g1_only <- setdiff(g1_hubs, g2_hubs)
        g2_only <- setdiff(g2_hubs, g1_hubs)
        
        jaccard_similarity <- length(shared_hubs) / length(union(g1_hubs, g2_hubs))
        
        new_row <- data.frame(
          Analysis_Type = "Genotype_Comparison",
          Condition1 = "G1",
          Condition2 = "G2",
          Group1_hubs = length(g1_hubs),
          Group2_hubs = length(g2_hubs),
          Shared_hubs = length(shared_hubs),
          Group1_unique = length(g1_only),
          Group2_unique = length(g2_only),
          Jaccard_similarity = jaccard_similarity,
          Comparison_Label = paste(tools::toTitleCase(tissue), "G1_vs_G2"),
          stringsAsFactors = FALSE
        )
        results <- rbind(results, new_row)
      }
    }
  }
  
  # Analyze overlap between tissues within genotypes
  for (genotype in c("g1", "g2")) {
    leaf_condition <- paste("leaf", genotype, sep = "_")
    root_condition <- paste("root", genotype, sep = "_")
    
    if (leaf_condition %in% names(hub_sets) && root_condition %in% names(hub_sets)) {
      leaf_hubs <- hub_sets[[leaf_condition]]
      root_hubs <- hub_sets[[root_condition]]
      
      if (length(leaf_hubs) > 0 && length(root_hubs) > 0) {
        shared_hubs <- intersect(leaf_hubs, root_hubs)
        leaf_only <- setdiff(leaf_hubs, root_hubs)
        root_only <- setdiff(root_hubs, leaf_hubs)
        
        jaccard_similarity <- length(shared_hubs) / length(union(leaf_hubs, root_hubs))
        
        new_row <- data.frame(
          Analysis_Type = "Tissue_Comparison",
          Condition1 = "Leaf",
          Condition2 = "Root",
          Group1_hubs = length(leaf_hubs),
          Group2_hubs = length(root_hubs),
          Shared_hubs = length(shared_hubs),
          Group1_unique = length(leaf_only),
          Group2_unique = length(root_only),
          Jaccard_similarity = jaccard_similarity,
          Comparison_Label = paste(toupper(genotype), "Leaf_vs_Root"),
          stringsAsFactors = FALSE
        )
        results <- rbind(results, new_row)
      }
    }
  }
  
  if (nrow(results) > 0) {
    cat("✓ Hub persistence analysis completed for", nrow(results), "comparisons\n")
    return(results)
  } else {
    cat("Warning: No valid hub comparisons could be performed\n")
    return(NULL)
  }
}

#' Track VIP Metabolites Over Time
track_vip_metabolites <- function(leaf_data, root_data, vip_data, shared_metabolites) {
  cat("Tracking VIP metabolites over time...\n")
  
  if (is.null(vip_data) || nrow(vip_data) == 0) {
    cat("Warning: No VIP data available\n")
    return(NULL)
  }
  
  # Get VIP metabolites that exist in our data
  vip_metabolites <- vip_data$Metabolite
  vip_in_data <- intersect(vip_metabolites, shared_metabolites)
  
  if (length(vip_in_data) == 0) {
    cat("Warning: No VIP metabolites found in shared metabolites\n")
    return(NULL)
  }
  
  cat("Tracking", length(vip_in_data), "VIP metabolites\n")
  
  # Combine leaf and root data
  # Find common columns between datasets
  common_cols <- intersect(names(leaf_data), names(root_data))
  required_cols <- c("Day", "Genotype", "Treatment")
  metabolite_cols <- intersect(shared_metabolites, common_cols)
  all_needed_cols <- unique(c(required_cols, metabolite_cols))
  
  # Combine data with only common columns
  combined_data <- rbind(
    leaf_data %>% select(all_of(all_needed_cols)) %>% mutate(Tissue = "Leaf"),
    root_data %>% select(all_of(all_needed_cols)) %>% mutate(Tissue = "Root")
  )
  
  # Create standardized result structure
  vip_tracking_results <- data.frame(
    Metabolite = character(),
    Genotype = character(),
    Tissue = character(),
    Treatment = numeric(),
    Temporal_slope = numeric(),
    Slope_p_value = numeric(),
    Trend_category = character(),
    CV = numeric(),
    Max_change = numeric(),
    N_timepoints = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Track each VIP metabolite
  for (metabolite in vip_in_data) {
    for (genotype in c("G1", "G2")) {
      for (tissue in c("Leaf", "Root")) {
        for (treatment in c(0, 1)) {
          subset_data <- combined_data %>%
            filter(Genotype == genotype, Tissue == tissue, Treatment == treatment) %>%
            select(Day, all_of(metabolite))
          
          if (nrow(subset_data) > 0 && metabolite %in% names(subset_data)) {
            # Calculate temporal trajectory
            temporal_summary <- subset_data %>%
              group_by(Day) %>%
              summarise(
                Mean_Value = mean(.data[[metabolite]], na.rm = TRUE),
                N_Observations = sum(!is.na(.data[[metabolite]])),
                .groups = 'drop'
              ) %>%
              filter(N_Observations > 0)
            
            if (nrow(temporal_summary) >= 2) {
              # Calculate temporal trend
              slope <- NA
              p_value <- NA
              
              if (nrow(temporal_summary) >= 3) {
                tryCatch({
                  lm_result <- lm(Mean_Value ~ Day, data = temporal_summary)
                  slope <- coef(lm_result)[2]
                  p_value <- summary(lm_result)$coefficients[2, 4]
                }, error = function(e) {
                  slope <<- NA
                  p_value <<- NA
                })
              } else {
                slope <- (temporal_summary$Mean_Value[2] - temporal_summary$Mean_Value[1]) / 
                  (temporal_summary$Day[2] - temporal_summary$Day[1])
              }
              
              # Classify trend
              trend_category <- case_when(
                is.na(slope) ~ "Unknown",
                abs(slope) < 0.01 ~ "Stable",
                slope > 0 ~ "Increasing",
                slope < 0 ~ "Decreasing",
                TRUE ~ "Stable"
              )
              
              # Calculate CV and max change
              cv_val <- NA
              max_change_val <- NA
              
              if (length(temporal_summary$Mean_Value) > 1) {
                cv_val <- sd(temporal_summary$Mean_Value) / mean(temporal_summary$Mean_Value)
                max_change_val <- max(temporal_summary$Mean_Value) - min(temporal_summary$Mean_Value)
              }
              
              # Create new row with consistent structure
              new_row <- data.frame(
                Metabolite = metabolite,
                Genotype = genotype,
                Tissue = tissue,
                Treatment = treatment,
                Temporal_slope = ifelse(is.na(slope), 0, slope),
                Slope_p_value = ifelse(is.na(p_value), 1, p_value),
                Trend_category = trend_category,
                CV = ifelse(is.na(cv_val), 0, cv_val),
                Max_change = ifelse(is.na(max_change_val), 0, max_change_val),
                N_timepoints = nrow(temporal_summary),
                stringsAsFactors = FALSE
              )
              
              # Add to results
              vip_tracking_results <- rbind(vip_tracking_results, new_row)
            }
          }
        }
      }
    }
  }
  
  if (nrow(vip_tracking_results) > 0) {
    cat("✓ VIP tracking completed for", nrow(vip_tracking_results), "conditions\n")
    return(vip_tracking_results)
  } else {
    cat("Warning: No VIP tracking results generated\n")
    return(NULL)
  }
}

#' Calculate Feature Change Counts
calculate_feature_changes <- function(leaf_data, root_data, shared_metabolites) {
  cat("Calculating feature change counts over time...\n")
  
  # Function to test temporal trends for each metabolite
  test_temporal_trends <- function(data, tissue_name) {
    results <- list()
    
    for (genotype in c("G1", "G2")) {
      for (treatment in c(0, 1)) {
        subset_data <- data %>%
          filter(Genotype == genotype, Treatment == treatment) %>%
          select(Day, all_of(shared_metabolites))
        
        if (nrow(subset_data) > 0) {
          # Test each metabolite for temporal trends
          metabolite_trends <- sapply(shared_metabolites, function(metabolite) {
            if (metabolite %in% names(subset_data)) {
              metabolite_data <- subset_data %>%
                select(Day, all_of(metabolite)) %>%
                filter(!is.na(.data[[metabolite]]))
              
              if (nrow(metabolite_data) >= 3) {
                # Perform linear regression
                lm_result <- tryCatch({
                  lm(get(metabolite) ~ Day, data = metabolite_data)
                }, error = function(e) NULL)
                
                if (!is.null(lm_result)) {
                  p_value <- summary(lm_result)$coefficients[2, 4]
                  return(p_value < 0.05)  # Significant trend
                }
              }
            }
            return(FALSE)
          })
          
          results[[paste(genotype, treatment, sep = "_")]] <- data.frame(
            Genotype = genotype,
            Treatment = treatment,
            Tissue = tissue_name,
            Significant_trends = sum(metabolite_trends, na.rm = TRUE),
            Total_features = length(shared_metabolites),
            Fraction_changing = sum(metabolite_trends, na.rm = TRUE) / length(shared_metabolites),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    return(do.call(rbind, results))
  }
  
  leaf_changes <- test_temporal_trends(leaf_data, "Leaf")
  root_changes <- test_temporal_trends(root_data, "Root")
  
  combined_changes <- rbind(leaf_changes, root_changes)
  cat("✓ Feature change analysis completed\n")
  return(combined_changes)
}

#' Calculate Feature Stability
calculate_feature_stability <- function(leaf_data, root_data, shared_metabolites) {
  cat("Calculating feature stability...\n")
  
  # Function to calculate resilient features at final timepoint
  calculate_resilient_features <- function(data, tissue_name) {
    results <- list()
    
    for (genotype in c("G1", "G2")) {
      # Compare final timepoint to initial (Treatment 1 vs 0)
      control_data <- data %>%
        filter(Genotype == genotype, Treatment == 0) %>%
        select(all_of(shared_metabolites)) %>%
        summarise_all(mean, na.rm = TRUE)
      
      stress_data <- data %>%
        filter(Genotype == genotype, Treatment == 1, Day == max(Day, na.rm = TRUE)) %>%
        select(all_of(shared_metabolites)) %>%
        summarise_all(mean, na.rm = TRUE)
      
      if (nrow(control_data) > 0 && nrow(stress_data) > 0) {
        # Calculate which features are "resilient" (little change)
        resilient_count <- 0
        total_count <- 0
        
        for (metabolite in shared_metabolites) {
          if (metabolite %in% names(control_data) && metabolite %in% names(stress_data)) {
            control_val <- control_data[[metabolite]]
            stress_val <- stress_data[[metabolite]]
            
            if (!is.na(control_val) && !is.na(stress_val) && control_val != 0) {
              fold_change <- abs(stress_val - control_val) / abs(control_val)
              if (fold_change < 0.2) {  # Less than 20% change considered resilient
                resilient_count <- resilient_count + 1
              }
              total_count <- total_count + 1
            }
          }
        }
        
        results[[genotype]] <- data.frame(
          Genotype = genotype,
          Tissue = tissue_name,
          Resilient_features = resilient_count,
          Total_features = total_count,
          Resilient_fraction = resilient_count / total_count * 100,
          stringsAsFactors = FALSE
        )
      }
    }
    
    return(do.call(rbind, results))
  }
  
  leaf_stability <- calculate_resilient_features(leaf_data, "Leaf")
  root_stability <- calculate_resilient_features(root_data, "Root")
  
  combined_stability <- rbind(leaf_stability, root_stability)
  cat("✓ Feature stability analysis completed\n")
  return(combined_stability)
}

#' Calculate Mean Response Validation Data
calculate_mean_responses <- function(leaf_data, root_data, shared_metabolites) {
  cat("Calculating mean response validation data...\n")
  
  process_tissue <- function(data, tissue_name) {
    data %>%
      filter(Treatment == 1) %>%  # Focus on stress treatment
      select(Genotype, all_of(shared_metabolites)) %>%
      pivot_longer(
        cols = all_of(shared_metabolites),
        names_to = "Metabolite",
        values_to = "Value"
      ) %>%
      group_by(Genotype) %>%
      summarise(
        Mean_Response = mean(abs(Value), na.rm = TRUE),
        Median_Response = median(abs(Value), na.rm = TRUE),
        SD_Response = sd(abs(Value), na.rm = TRUE),
        N_Observations = sum(!is.na(Value)),
        .groups = 'drop'
      ) %>%
      mutate(
        Tissue = tissue_name,
        SE_Response = SD_Response / sqrt(N_Observations),
        CI95_Lower = Mean_Response - 1.96 * SE_Response,
        CI95_Upper = Mean_Response + 1.96 * SE_Response
      )
  }
  
  leaf_responses <- process_tissue(leaf_data, "Leaf")
  root_responses <- process_tissue(root_data, "Root")
  
  combined_responses <- bind_rows(leaf_responses, root_responses)
  cat("✓ Calculated mean responses for", nrow(combined_responses), "conditions\n")
  return(combined_responses)
}

#' Generate Null Distribution for Validation
generate_null_distributions <- function(leaf_data, root_data, shared_metabolites, n_permutations = 500) {
  cat("Generating null distributions (", n_permutations, "permutations)...\n")
  
  null_results <- list()
  
  # Sample a subset of metabolites for computational efficiency
  sample_metabolites <- sample(shared_metabolites, min(50, length(shared_metabolites)))
  
  for (i in 1:n_permutations) {
    if (i %% 100 == 0) cat("  Permutation", i, "/", n_permutations, "\n")
    
    for (genotype in c("G1", "G2")) {
      for (tissue_data in list(list(data = leaf_data, name = "Leaf"), 
                               list(data = root_data, name = "Root"))) {
        
        tissue_name <- tissue_data$name
        data <- tissue_data$data
        
        # Get data for this genotype
        genotype_data <- data %>%
          filter(Genotype == genotype) %>%
          select(Treatment, all_of(sample_metabolites))
        
        if (nrow(genotype_data) > 10) {
          # Randomly permute treatment labels
          permuted_data <- genotype_data
          permuted_data$Treatment <- sample(genotype_data$Treatment)
          
          # Calculate "effect size" for one representative metabolite
          metabolite <- sample_metabolites[1]
          control_vals <- permuted_data[permuted_data$Treatment == 0, metabolite]
          treatment_vals <- permuted_data[permuted_data$Treatment == 1, metabolite]
          
          control_clean <- control_vals[!is.na(control_vals)]
          treatment_clean <- treatment_vals[!is.na(treatment_vals)]
          
          if (length(control_clean) >= 3 && length(treatment_clean) >= 3) {
            tryCatch({
              cliff_result <- cliff.delta(treatment_clean, control_clean)
              
              null_results[[paste(i, genotype, tissue_name, sep = "_")]] <- data.frame(
                Permutation = i,
                Genotype = genotype,
                Tissue = tissue_name,
                Null_Effect_Size = cliff_result$estimate
              )
            }, error = function(e) {
              # Skip failed permutations
            })
          }
        }
      }
    }
  }
  
  null_df <- do.call(rbind, null_results)
  cat("✓ Generated", nrow(null_df), "null effect size estimates\n")
  return(null_df)
}

# 7. EXECUTE ANALYSES
# =============================================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("EXECUTING COMPREHENSIVE NETWORK & MECHANISM ANALYSES\n")
cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")

# Run all analyses
effect_sizes_data <- calculate_effect_sizes(leaf_data, root_data, shared_metabolites)
hub_persistence_data <- calculate_hub_persistence(hub_data)
vip_tracking_data <- track_vip_metabolites(leaf_data, root_data, vip_data, shared_metabolites)
feature_changes_data <- calculate_feature_changes(leaf_data, root_data, shared_metabolites)
feature_stability_data <- calculate_feature_stability(leaf_data, root_data, shared_metabolites)
mean_responses_data <- calculate_mean_responses(leaf_data, root_data, shared_metabolites)
null_distributions_data <- generate_null_distributions(leaf_data, root_data, shared_metabolites, 500)

# 8. PANEL GENERATION FUNCTIONS
# =============================================================================

#' Panel A: Hub Persistence Analysis
create_panel_a <- function(hub_persistence_data, plot_sizing, color_palette) {
  cat("Creating Panel A: Hub persistence analysis...\n")
  
  if (is.null(hub_persistence_data) || nrow(hub_persistence_data) == 0) {
    # Create placeholder plot if no data
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Hub data not available", size = plot_sizing$placeholder_text_size) +
      theme_void() +
      labs(title = "Hub Persistence Analysis", tag = "A")
    return(p)
  }
  
  # Create Jaccard similarity plot
  p <- ggplot(hub_persistence_data, aes(x = reorder(Comparison_Label, Jaccard_similarity), 
                                        y = Jaccard_similarity)) +
    geom_col(aes(fill = Analysis_Type), width = 0.7, alpha = 0.8) +
    geom_text(aes(label = round(Jaccard_similarity, 3)), 
              hjust = -0.1, fontface = "plain", size = plot_sizing$geom_text_size) +
    
    scale_fill_manual(values = c("Genotype_Comparison" = color_palette$genotype[["G1"]], 
                                 "Tissue_Comparison" = color_palette$genotype[["G2"]]), 
                      name = "Analysis Type",
                      labels = c("Genotype_Comparison" = "Genotype", 
                                 "Tissue_Comparison" = "Tissue")) +
    
    scale_y_continuous(limits = c(0, max(hub_persistence_data$Jaccard_similarity) * 1.3),
                       labels = number_format(accuracy = 0.001)) +
    
    coord_flip() +
    
    labs(
      title = "Hub Persistence Analysis",
      subtitle = "Jaccard similarity between conditions",
      x = "Comparison",
      y = "Jaccard Similarity",
      tag = "A"
    ) +
    
    create_nature_theme(plot_sizing, y_title_margin = 25) +  # Larger margin for Panel A
    theme(axis.text.y = element_text(size = plot_sizing$axis_y_text_size))
  
  return(p)
}

#' Panel B: VIP Metabolite Tracking
create_panel_b <- function(vip_tracking_data, plot_sizing, color_palette) {
  cat("Creating Panel B: VIP metabolite tracking...\n")
  
  if (is.null(vip_tracking_data)) {
    # Create placeholder plot if no data
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "VIP tracking data not available", size = plot_sizing$placeholder_text_size) +
      theme_void() +
      labs(title = "VIP Metabolite Tracking", tag = "B")
    return(p)
  }
  
  # Focus on stress treatment
  stress_vip <- vip_tracking_data %>%
    filter(Treatment == 1) %>%
    filter(!is.na(Temporal_slope))
  
  if (nrow(stress_vip) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No stress VIP data available", size = plot_sizing$placeholder_text_size) +
      theme_void() +
      labs(title = "VIP Metabolite Tracking", tag = "B")
    return(p)
  }
  
  p <- ggplot(stress_vip, aes(x = Tissue, y = Temporal_slope, fill = Genotype, color = Genotype)) +
    geom_violin(alpha = 0.7, position = position_dodge(0.8)) +
    geom_boxplot(width = 0.2, position = position_dodge(0.8), 
                 outlier.shape = NA, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    
    scale_fill_manual(values = color_palette$genotype, name = "Genotype", labels = c("G1" = "Genotype 1", "G2" = "Genotype 2")) +
    scale_color_manual(values = color_palette$genotype, name = "Genotype", labels = c("G1" = "Genotype 1", "G2" = "Genotype 2")) +
    
    labs(
      title = "VIP Metabolite Tracking",
      subtitle = "Temporal slope distribution under stress",
      x = "Tissue",
      y = "Temporal Slope (Change/Day)",
      tag = "B"
    ) +
    
    create_nature_theme(plot_sizing, y_title_margin = 8)  # Smaller margin for Panel B
  
  return(p)
}

#' Panel C: Feature Change Counts
create_panel_c <- function(feature_changes_data, plot_sizing, color_palette) {
  cat("Creating Panel C: Feature change counts...\n")
  
  if (is.null(feature_changes_data)) {
    # Create placeholder plot if no data
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Feature change data not available", size = plot_sizing$placeholder_text_size) +
      theme_void() +
      labs(title = "Feature Change Counts", tag = "C")
    return(p)
  }
  
  # Focus on stress treatment
  stress_changes <- feature_changes_data %>%
    filter(Treatment == 1)
  
  p <- ggplot(stress_changes, aes(x = Tissue, y = Significant_trends, fill = Genotype)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(label = Significant_trends), 
              position = position_dodge(0.8), vjust = -0.2,  # Changed from -0.5 to -0.2
              fontface = "plain", size = plot_sizing$geom_text_size) +
    
    scale_fill_manual(values = color_palette$genotype, name = "Genotype") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Added padding at top
    
    labs(
      title = "Feature Change Counts",
      subtitle = "Number of features with significant temporal trends",
      x = "Tissue",
      y = "Features with Significant Trends",
      tag = "C"
    ) +
    
    create_nature_theme(plot_sizing, y_title_margin = 8)  # Smaller margin for Panel C
  
  return(p)
}

#' Panel D: Feature Stability
create_panel_d <- function(feature_stability_data, plot_sizing, color_palette) {
  cat("Creating Panel D: Feature stability...\n")
  
  if (is.null(feature_stability_data)) {
    # Create placeholder plot if no data
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Feature stability data not available", size = plot_sizing$placeholder_text_size) +
      theme_void() +
      labs(title = "Feature Stability", tag = "D")
    return(p)
  }
  
  p <- ggplot(feature_stability_data, aes(x = Tissue, y = Resilient_fraction, fill = Genotype)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(label = paste0(round(Resilient_fraction, 1), "%")), 
              position = position_dodge(0.8), vjust = -0.2,  # Changed from -0.5 to -0.2
              fontface = "plain", size = plot_sizing$geom_text_size) +
    
    scale_fill_manual(values = color_palette$genotype, name = "Genotype") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.1))) +  # Added padding at top
    
    labs(
      title = "Feature Stability",
      subtitle = "Resilient features at final timepoint",
      x = "Tissue",
      y = "Resilient Features (%)",
      tag = "D"
    ) +
    
    create_nature_theme(plot_sizing, y_title_margin = 8)  # Smaller margin for Panel D
  
  return(p)
}

#' Panel E: Mean Response Validation
create_panel_e <- function(mean_responses_data, plot_sizing, color_palette) {
  cat("Creating Panel E: Mean response validation...\n")
  
  p <- ggplot(mean_responses_data, aes(x = Genotype, y = Mean_Response, fill = Tissue)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, CI95_Lower), ymax = CI95_Upper), 
                  position = position_dodge(0.8), width = 0.25) +
    
    # Add value labels on bars
    geom_text(aes(label = round(Mean_Response, 0), y = CI95_Upper), 
              position = position_dodge(0.8), vjust = -0.2, size = plot_sizing$geom_text_size) +  # Changed from -0.5 to -0.2
    
    # Styling with original colors
    scale_fill_manual(values = color_palette$tissue, name = "Tissue") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Added padding at top
    
    labs(
      title = "Mean Response Validation",
      subtitle = "Static comparison of tissue-specific response intensities",
      x = "Genotype",
      y = "Mean Molecular Response (Arbitrary Units)",
      tag = "E"
    ) +
    
    create_nature_theme(plot_sizing, y_title_margin = 8)  # Smaller margin for Panel E
  
  return(p)
}

#' Panel F: Null Distribution Control
create_panel_f <- function(null_distributions_data, effect_sizes_data, plot_sizing, color_palette) {
  cat("Creating Panel F: Null distribution validation...\n")
  
  # Calculate observed effect sizes summary
  observed_summary <- effect_sizes_data %>%
    group_by(Genotype, Tissue) %>%
    summarise(
      Mean_Effect = mean(Effect_Size, na.rm = TRUE),
      SD_Effect = sd(Effect_Size, na.rm = TRUE),
      Type = "Observed",
      .groups = 'drop'
    )
  
  # Calculate null distribution summary
  null_summary <- null_distributions_data %>%
    group_by(Genotype, Tissue) %>%
    summarise(
      Mean_Effect = mean(Null_Effect_Size, na.rm = TRUE),
      SD_Effect = sd(Null_Effect_Size, na.rm = TRUE),
      Type = "Null",
      .groups = 'drop'
    )
  
  combined_data <- bind_rows(observed_summary, null_summary) %>%
    mutate(Comparison = paste(Genotype, Tissue, sep = "_"))
  
  p <- ggplot(combined_data, aes(x = Comparison, y = abs(Mean_Effect), fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs(Mean_Effect) - SD_Effect), 
                      ymax = abs(Mean_Effect) + SD_Effect), 
                  position = position_dodge(0.8), width = 0.25) +
    
    # Styling
    scale_fill_manual(values = c("Observed" = color_palette$accent, 
                                 "Null" = color_palette$neutral[2]), 
                      name = "Distribution") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Ensure y-axis starts at 0
    
    labs(
      title = "Null Distribution Validation",
      subtitle = "Methodological control: observed vs. permuted effect sizes",
      x = "Genotype × Tissue",
      y = "Absolute Mean Effect Size",
      tag = "F"
    ) +
    
    create_nature_theme(plot_sizing, y_title_margin = 8) +  # Smaller margin for Panel F
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# 9. GENERATE PANELS
# =============================================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("GENERATING SIX-PANEL FIGURE (A-F)\n")
cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")

panel_a <- create_panel_a(hub_persistence_data, plot_sizing, color_palette)
panel_b <- create_panel_b(vip_tracking_data, plot_sizing, color_palette)
panel_c <- create_panel_c(feature_changes_data, plot_sizing, color_palette)
panel_d <- create_panel_d(feature_stability_data, plot_sizing, color_palette)
panel_e <- create_panel_e(mean_responses_data, plot_sizing, color_palette)
panel_f <- create_panel_f(null_distributions_data, effect_sizes_data, plot_sizing, color_palette)

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
safe_ggsave_png(panel_a, "Fig2_Panel_A_HubPersistence", plot_sizing$individual_panel_width, plot_sizing$individual_panel_height)
safe_ggsave_png(panel_b, "Fig2_Panel_B_VIPTracking", plot_sizing$individual_panel_width, plot_sizing$individual_panel_height)
safe_ggsave_png(panel_c, "Fig2_Panel_C_FeatureChanges", plot_sizing$individual_panel_width, plot_sizing$individual_panel_height)
safe_ggsave_png(panel_d, "Fig2_Panel_D_FeatureStability", plot_sizing$individual_panel_width, plot_sizing$individual_panel_height)
safe_ggsave_png(panel_e, "Fig2_Panel_E_MeanResponseValidation", plot_sizing$individual_panel_width, plot_sizing$individual_panel_height)
safe_ggsave_png(panel_f, "Fig2_Panel_F_NullDistribution", plot_sizing$individual_panel_width, plot_sizing$individual_panel_height)

# 11. CREATE COMBINED SIX-PANEL FIGURE
# =============================================================================
cat("Creating combined six-panel figure...\n")

# Combine panels in 2x3 grid with reduced spacing
combined_figure <- plot_grid(
  panel_a, panel_b, panel_c,
  panel_d, panel_e, panel_f,
  ncol = 3, nrow = 2,
  align = 'hv',
  hjust = 0, vjust = 1,  # Alignment settings
  rel_widths = c(1, 1, 1),  # Equal widths
  rel_heights = c(1, 1)     # Equal heights
)

# Save combined figure
safe_ggsave_png(combined_figure, "Fig2_SixPanel_AF_Combined", plot_sizing$combined_panel_width, plot_sizing$combined_panel_height)

# 12. SAVE PROCESSED DATA
# =============================================================================
cat("Saving processed data...\n")

# Save data if it exists
if (!is.null(effect_sizes_data)) {
  write.csv(effect_sizes_data, file.path(output_dir, "processed_effect_sizes_sixpanel.csv"), row.names = FALSE)
}
if (!is.null(hub_persistence_data)) {
  write.csv(hub_persistence_data, file.path(output_dir, "processed_hub_persistence_data_sixpanel.csv"), row.names = FALSE)
}
if (!is.null(vip_tracking_data)) {
  write.csv(vip_tracking_data, file.path(output_dir, "processed_vip_tracking_data_sixpanel.csv"), row.names = FALSE)
}
if (!is.null(feature_changes_data)) {
  write.csv(feature_changes_data, file.path(output_dir, "processed_feature_changes_data_sixpanel.csv"), row.names = FALSE)
}
if (!is.null(feature_stability_data)) {
  write.csv(feature_stability_data, file.path(output_dir, "processed_feature_stability_data_sixpanel.csv"), row.names = FALSE)
}
if (!is.null(mean_responses_data)) {
  write.csv(mean_responses_data, file.path(output_dir, "processed_mean_responses_sixpanel.csv"), row.names = FALSE)
}
if (!is.null(null_distributions_data)) {
  write.csv(null_distributions_data, file.path(output_dir, "processed_null_distributions_sixpanel.csv"), row.names = FALSE)
}

# 13. SUMMARY STATISTICS
# =============================================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
cat("SIX-PANEL FIGURE ANALYSIS SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")

cat("Analysis Results:\n")
cat("• Effect sizes calculated:", ifelse(is.null(effect_sizes_data), "No data", nrow(effect_sizes_data)), "metabolite-tissue-genotype combinations\n")
cat("• Hub persistence analysis:", ifelse(is.null(hub_persistence_data), "No data", nrow(hub_persistence_data)), "comparisons\n")
cat("• VIP tracking analysis:", ifelse(is.null(vip_tracking_data), "No data", nrow(vip_tracking_data)), "conditions\n")
cat("• Feature changes analysis:", ifelse(is.null(feature_changes_data), "No data", nrow(feature_changes_data)), "conditions\n")
cat("• Feature stability analysis:", ifelse(is.null(feature_stability_data), "No data", nrow(feature_stability_data)), "conditions\n")
cat("• Mean response analysis:", ifelse(is.null(mean_responses_data), "No data", nrow(mean_responses_data)), "conditions\n")
cat("• Null distribution permutations:", ifelse(is.null(null_distributions_data), "No data", max(null_distributions_data$Permutation, na.rm = TRUE)), "\n")
cat("• Shared metabolites used:", length(shared_metabolites), "\n")
cat("• Output files generated: 14\n")

cat("\n✅ Six-Panel Figure Script completed successfully!\n")
cat("📁 All outputs saved to:", output_dir, "\n")
cat("🎯 Combined figure shows all six panels (A-F) in a 2×3 grid layout\n")