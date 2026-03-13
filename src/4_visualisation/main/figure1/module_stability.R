#-------------------------------------------------------------------------------
# SCRIPT: module_stability.R
#
# DESCRIPTION:
# This script assesses the stability of metabolic modules within different
# plant tissues (leaf and root) and genotypes (G1, G2). It uses a
# permutation-based approach to calculate two key metrics:
#   1. Module Preservation: The mean absolute correlation within a module.
#   2. Module Coherence: The proportion of correlations above a defined
#      threshold (e.g., 0.7).
#
# The script generates a scatter plot visualizing these metrics for each
# tissue-genotype combination and saves the plot and underlying data.
#
# INPUTS:
#   - Leaf data CSV: Contains metabolic features for leaf samples.
#   - Root data CSV: Contains metabolic features for root samples.
#   - VIP scores CSV: Used for feature selection (though not directly in this
#     script, its path is defined).
#
# OUTPUTS:
#   - "module_stability.pdf": A PDF scatter plot of module stability metrics.
#   - "module_stability.png": A PNG scatter plot of module stability metrics.
#   - "module_stability_metrics.csv": A CSV file with the calculated stability
#     metrics.
#
#
#
# USAGE:
# Ensure the file paths in the 'CONFIGURATION' section are correct, then
# execute the script in an R environment with ggplot2 installed.
#
#-------------------------------------------------------------------------------

#===============================================================================
# 1. SETUP
#===============================================================================

# Load necessary libraries
# Ensure ggplot2 is installed: install.packages("ggplot2")
library(ggplot2)

# Set seed for reproducibility across the entire script
set.seed(42)


#===============================================================================
# 2. CONFIGURATION
#===============================================================================
# Update these paths before running the script.

# Input file paths
leaf_data_path <- "C:/Users/USER/Desktop/data_chem_3_10/data/data/n_p_l.csv"
root_data_path <- "C:/Users/USER/Desktop/data_chem_3_10/data/data/n_p_r.csv"
vip_path       <- "C:/Users/USER/Desktop/data_chem_3_10/output/results/vip_bonferroni/VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv"

# Output directory
output_dir <- "C:/Users/USER/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000/module_stability"

# Analysis parameters
n_permutations      <- 1000 # Number of permutations for the significance test
coherence_threshold <- 0.7  # Correlation threshold for coherence calculation


#===============================================================================
# 3. DATA LOADING AND PREPARATION
#===============================================================================

# Create the output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load the datasets
# Using tryCatch for more robust error handling during file loading
tryCatch({
    leaf_data <- read.csv(leaf_data_path)
    root_data <- read.csv(root_data_path)
    vip_data  <- read.csv(vip_path) # This is loaded but not used in this script
}, error = function(e) {
    stop("Error loading data files. Please check the paths in the CONFIGURATION section. Details: ", e$message)
})


#===============================================================================
# 4. FUNCTION DEFINITION
#===============================================================================

#' Calculate Module Stability Metrics
#'
#' This function computes module preservation and coherence for a given dataset
#' and genotype, and assesses their statistical significance via permutation testing.
#'
#' @param data A data frame containing metabolic features and a 'Genotype' column.
#' @param genotype A character string specifying the genotype to analyze (e.g., "G1").
#' @param n_perms An integer for the number of permutations to perform.
#' @param coherence_thresh A numeric value for the correlation coherence threshold.
#'
#' @return A numeric vector containing:
#'   1. Observed Module Preservation
#'   2. Observed Module Coherence
#'   3. p-value for Preservation
#'   4. p-value for Coherence
#'
calculate_module_metrics <- function(data, genotype, n_perms, coherence_thresh) {
  # Set seed within the function for consistent permutation results
  set.seed(42)

  # Select feature columns (starting with N_Cluster_ or P_Cluster_)
  feature_cols <- grep("^(N_Cluster_|P_Cluster_)", names(data), value = TRUE)
  genotype_data <- data[data$Genotype == genotype, feature_cols]

  # --- Calculate observed metrics ---
  correlation_matrix <- cor(genotype_data, use = "pairwise.complete.obs")

  # Module Preservation: Mean absolute correlation
  observed_preservation <- mean(abs(correlation_matrix), na.rm = TRUE)

  # Module Coherence: Proportion of correlations > threshold
  observed_coherence <- mean(correlation_matrix > coherence_thresh, na.rm = TRUE)

  # --- Perform permutation test to establish a null distribution ---
  null_preservation <- numeric(n_perms)
  null_coherence    <- numeric(n_perms)

  for (i in 1:n_perms) {
    # Permute each feature column independently to break correlation structures
    permuted_data <- apply(genotype_data, 2, sample)
    perm_cor      <- cor(permuted_data, use = "pairwise.complete.obs")

    null_preservation[i] <- mean(abs(perm_cor), na.rm = TRUE)
    null_coherence[i]    <- mean(perm_cor > coherence_thresh, na.rm = TRUE)
  }

  # Calculate p-values: the proportion of null values as or more extreme than observed
  p_preservation <- mean(null_preservation >= observed_preservation)
  p_coherence    <- mean(null_coherence >= observed_coherence)

  return(c(
    Module_Preservation = observed_preservation,
    Module_Coherence    = observed_coherence,
    P_Preservation      = p_preservation,
    P_Coherence         = p_coherence
  ))
}


#===============================================================================
# 5. ANALYSIS
#===============================================================================

# Initialize a list to store results from each analysis run
results_list <- list()

# Define genotypes to iterate over
genotypes <- c("G1", "G2")

# Process each tissue and genotype
for (geno in genotypes) {
  # Analyze leaf data
  leaf_metrics <- calculate_module_metrics(leaf_data, geno, n_permutations, coherence_threshold)
  results_list[[length(results_list) + 1]] <- data.frame(
    Tissue = "Leaf",
    Genotype = geno,
    as.list(leaf_metrics)
  )

  # Analyze root data
  root_metrics <- calculate_module_metrics(root_data, geno, n_permutations, coherence_threshold)
  results_list[[length(results_list) + 1]] <- data.frame(
    Tissue = "Root",
    Genotype = geno,
    as.list(root_metrics)
  )
}

# Combine all results into a single data frame
stability_metrics <- do.call(rbind, results_list)

# Create a combined 'Group' column for plotting
stability_metrics$Group <- paste(stability_metrics$Tissue, stability_metrics$Genotype)


#===============================================================================
# 6. VISUALIZATION
#===============================================================================

# Define a color palette for the plot for consistency
color_palette <- c("Leaf G1" = "#41e085",
                   "Leaf G2" = "#209150",
                   "Root G1" = "#33d6d3",
                   "Root G2" = "#1e597d")

# Create the scatter plot using ggplot2
stability_plot <- ggplot(stability_metrics,
                         aes(x = Module_Preservation, y = Module_Coherence, color = Group)) +
  geom_point(size = 6, alpha = 0.7) +
  scale_color_manual(name = "Tissue-Genotype", values = color_palette) +
  labs(
    title = "Module Stability Analysis",
    subtitle = paste("Based on", formatC(n_permutations, big.mark=","), "permutations"),
    x = "Module Preservation Score",
    y = "Module Coherence Score"
  ) +
  theme_minimal(base_size = 16) + # Set a base font size for the theme
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA)
  )


#===============================================================================
# 7. SAVE OUTPUTS
#===============================================================================

# Construct full file paths
output_plot_pdf <- file.path(output_dir, "module_stability.pdf")
output_plot_png <- file.path(output_dir, "module_stability.png")
output_data_csv <- file.path(output_dir, "module_stability_metrics.csv")

# Save the plot in both PDF and PNG formats
ggsave(output_plot_pdf, stability_plot, width = 8, height = 6, device = "pdf")
ggsave(output_plot_png, stability_plot, width = 8, height = 6, dpi = 300)

# Save the metrics data to a CSV file
write.csv(stability_metrics, output_data_csv, row.names = FALSE)

# Print a confirmation message to the console
cat("Script finished.\n")
cat("Outputs saved to:", output_dir, "\n")
