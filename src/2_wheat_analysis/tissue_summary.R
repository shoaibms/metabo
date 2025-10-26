########################################
#####summarise fdr_results.csv #########
########################################


library(tidyverse)
library(data.table)

# Read FDR results
fdr_results <- fread(file.path(data_dir, "fdr_results.csv"))

# Create summary table
summary_table <- fdr_results %>%
  # Filter for significant results
  filter(fdr < 0.05) %>%
  # Group by tissue type
  group_by(Tissue.type) %>%
  # Create summary
  mutate(
    sig_level = case_when(
      fdr < 0.001 ~ "***",
      fdr < 0.01 ~ "**",
      fdr < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(Tissue.type, fdr) %>%
  select(
    Metabolite,
    Tissue = Tissue.type,
    `P_value` = p_value,
    FDR = fdr,
    Significance = sig_level
  )

# Create statistical overview
stat_overview <- fdr_results %>%
  group_by(Tissue.type) %>%
  summarise(
    total_metabolites = n(),
    significant = sum(fdr < 0.05),
    highly_significant = sum(fdr < 0.01),
    percent_significant = round(100 * sum(fdr < 0.05) / n(), 1)
  )

# Save results as CSV files
write.csv(summary_table, 
          file = file.path(save_dir, "FDR_Significant_Results.csv"), 
          row.names = FALSE)

write.csv(stat_overview, 
          file = file.path(save_dir, "FDR_Statistical_Overview.csv"), 
          row.names = FALSE)

# Print overview to console
print("Statistical Overview:")
print(stat_overview)

# Return the summary table invisibly
invisible(summary_table)



########################################
#####summarise  #######################
########################################

library(dplyr)
library(tidyr)
library(data.table)

# Read effect sizes data
data_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
save_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"

effect_data <- fread(file.path(data_dir, "effect_sizes.csv"))

# Create summary table
summary_table <- effect_data %>%
  group_by(Tissue.type, Genotype) %>%
  summarize(
    n_metabolites = n(),
    mean_effect = round(mean(effect_size, na.rm = TRUE), 3),
    median_effect = round(median(effect_size, na.rm = TRUE), 3),
    ci_lower = round(mean(ci_lower, na.rm = TRUE), 3),
    ci_upper = round(mean(ci_upper, na.rm = TRUE), 3),
    large_effects = sum(magnitude == "large", na.rm = TRUE),
    medium_effects = sum(magnitude == "medium", na.rm = TRUE),
    small_effects = sum(magnitude == "small", na.rm = TRUE),
    negligible_effects = sum(magnitude == "negligible", na.rm = TRUE),
    pct_significant = round(sum(abs(effect_size) > 0.3, na.rm = TRUE) / n() * 100, 1)
  ) %>%
  mutate(
    Tissue = ifelse(Tissue.type == "L", "Leaf", "Root"),
    effect_distribution = sprintf("%d/%d/%d/%d", 
                                  large_effects, medium_effects, 
                                  small_effects, negligible_effects)
  ) %>%
  select(Tissue, Genotype, n_metabolites, mean_effect, median_effect, 
         ci_lower, ci_upper, effect_distribution, pct_significant)

# Format as publication-ready table
formatted_table <- data.frame(
  "Tissue" = summary_table$Tissue,
  "Genotype" = summary_table$Genotype,
  "N" = summary_table$n_metabolites,
  "Mean Effect (95% CI)" = sprintf("%.3f (%.3f to %.3f)", 
                                   summary_table$mean_effect,
                                   summary_table$ci_lower,
                                   summary_table$ci_upper),
  "Median Effect" = summary_table$median_effect,
  "Effect Distribution*" = summary_table$effect_distribution,
  "% Significant†" = summary_table$pct_significant
)

# Save table
write.csv(formatted_table, 
          file.path(save_dir, "Table_effect_sizes_summary.csv"), 
          row.names = FALSE)

# Print table structure for manuscript
cat("Table 1: Tissue-specific metabolic effect sizes under osmotic stress\n\n")
print(formatted_table)
cat("\n*Effect Distribution: Format is Large/Medium/Small/Negligible effects")
cat("\n†Percentage of metabolites with |effect size| > 0.3")






########################################
#####summarise effect_size.csv #########
########################################

library(dplyr)
library(data.table)

# Read data
data_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
save_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
effect_data <- fread(file.path(data_dir, "effect_sizes.csv"))

# Create improved summary with genotype
summary_stats <- effect_data %>%
  group_by(Tissue.type, Genotype) %>%
  summarize(
    total_metabolites = n(),
    # Effect sizes
    mean_effect = round(mean(abs(effect_size), na.rm = TRUE), 3),
    se_effect = round(sd(abs(effect_size), na.rm = TRUE)/sqrt(n()), 3),
    # Response distribution
    large_effects = sum(magnitude == "large", na.rm = TRUE),
    medium_effects = sum(magnitude == "medium", na.rm = TRUE),
    small_effects = sum(magnitude == "small", na.rm = TRUE),
    # Percentages
    pct_large = round(large_effects/n() * 100, 1),
    pct_medium = round(medium_effects/n() * 100, 1),
    pct_small = round(small_effects/n() * 100, 1)
  ) %>%
  mutate(
    Tissue = ifelse(Tissue.type == "L", "Leaf", "Root"),
    effect_distribution = sprintf("%d/%d/%d", 
                                  large_effects, medium_effects, 
                                  small_effects)
  ) %>%
  select(Tissue, Genotype, total_metabolites, mean_effect, se_effect,
         effect_distribution, pct_large, pct_medium, pct_small)

# Also create tissue-only summary
tissue_summary <- effect_data %>%
  group_by(Tissue.type) %>%
  summarize(
    total_metabolites = n(),
    mean_effect = round(mean(abs(effect_size), na.rm = TRUE), 3),
    se_effect = round(sd(abs(effect_size), na.rm = TRUE)/sqrt(n()), 3),
    large_effects = sum(magnitude == "large", na.rm = TRUE),
    medium_effects = sum(magnitude == "medium", na.rm = TRUE),
    small_effects = sum(magnitude == "small", na.rm = TRUE),
    pct_large = round(large_effects/n() * 100, 1),
    pct_medium = round(medium_effects/n() * 100, 1),
    pct_small = round(small_effects/n() * 100, 1)
  ) %>%
  mutate(
    Tissue = ifelse(Tissue.type == "L", "Leaf", "Root"),
    effect_distribution = sprintf("%d/%d/%d", 
                                  large_effects, medium_effects, 
                                  small_effects)
  ) %>%
  select(Tissue, total_metabolites, mean_effect, se_effect,
         effect_distribution, pct_large, pct_medium, pct_small)

# Save both summaries
write.csv(summary_stats, 
          file.path(save_dir, "effect_sizes_summary_by_genotype.csv"), 
          row.names = FALSE)
write.csv(tissue_summary, 
          file.path(save_dir, "effect_sizes_summary_by_tissue.csv"), 
          row.names = FALSE)

# Print both summaries
cat("\nSummary by Tissue and Genotype:\n")
print(summary_stats)
cat("\nSummary by Tissue:\n")
print(tissue_summary)
# Save formatted table
write.csv(summary_stats, 
          file.path(save_dir, "Table_effect_sizes_tissue_summary.csv"), 
          row.names = FALSE)

# Print formatted table
print(summary_stats)




#########################################
##### summarise temporal_corr .csv ######
#########################################

library(dplyr)
library(tidyr)
library(data.table)

# Set paths
data_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
save_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"

# Read data
corr_data <- fread(file.path(data_dir, "temporal_corr.csv"))

# Function to process correlations
process_correlations <- function(df, corr_type) {
  if(corr_type == "spearman") {
    p_col <- "spearman_p"
    corr_col <- "spearman_rho"
  } else {
    p_col <- "kendall_p"
    corr_col <- "kendall_tau"
  }
  
  total_n <- nrow(df)
  
  # Filter significant correlations
  sig_data <- df[df[[p_col]] < 0.05, ]
  
  # Calculate percentages
  sig_pct <- nrow(sig_data) / total_n * 100
  pos_pct <- sum(sig_data[[corr_col]] > 0) / total_n * 100
  neg_pct <- sum(sig_data[[corr_col]] < 0) / total_n * 100
  
  # Get ranges for positive and negative correlations
  pos_range <- if(any(sig_data[[corr_col]] > 0)) {
    pos_vals <- sig_data[[corr_col]][sig_data[[corr_col]] > 0]
    sprintf("%.2f to %.2f", min(pos_vals), max(pos_vals))
  } else {
    "NA"
  }
  
  neg_range <- if(any(sig_data[[corr_col]] < 0)) {
    neg_vals <- sig_data[[corr_col]][sig_data[[corr_col]] < 0]
    sprintf("%.2f to %.2f", min(neg_vals), max(neg_vals))
  } else {
    "NA"
  }
  
  data.frame(
    Significant_pct = round(sig_pct, 1),
    Positive_pct = round(pos_pct, 1),
    Negative_pct = round(neg_pct, 1),
    Positive_range = pos_range,
    Negative_range = neg_range
  )
}

# Generate summaries
spearman_summary <- corr_data %>%
  group_by(Tissue.type, Genotype) %>%
  group_modify(~process_correlations(., "spearman")) %>%
  ungroup()

kendall_summary <- corr_data %>%
  group_by(Tissue.type, Genotype) %>%
  group_modify(~process_correlations(., "kendall")) %>%
  ungroup()

# Format tables
format_table <- function(data) {
  data %>%
    rename(
      "Tissue" = Tissue.type,
      "Significant (%)" = Significant_pct,
      "Positive (%)" = Positive_pct,
      "Negative (%)" = Negative_pct,
      "Positive Range" = Positive_range,
      "Negative Range" = Negative_range
    )
}

spearman_table <- format_table(spearman_summary)
kendall_table <- format_table(kendall_summary)

# Save tables
write.csv(spearman_table, file.path(save_dir, "spearman_correlation_table.csv"), row.names = FALSE)
write.csv(kendall_table, file.path(save_dir, "kendall_correlation_table.csv"), row.names = FALSE)

# Print tables
print("Spearman Correlation Summary:")
print(spearman_table, row.names = FALSE)
print("\nKendall Correlation Summary:")