# SCRIPT: tissue_plot.R
# AUTHOR: [Your Name/Lab Name]
# DATE: [Date of last revision]
#
# DESCRIPTION:
# This script generates a two-panel plot to visualize the effect sizes of
# metabolic changes between two genotypes (G1 and G2) in different tissue types
# (Leaf and Root).
#
# Panel 1: A ridge plot showing the distribution of effect sizes (Cliff's Delta)
#          for each genotype and tissue type.
# Panel 2: A bar chart comparing the ratio of significant metabolic changes
#          (Leaf vs. Root) for different effect size categories.
#
# The final combined plot is saved as both a PDF and a PNG file.

# --- 1. Load Libraries ---
# Using explicit `library()` calls for clarity.
library(ggplot2)    # For creating plots
library(dplyr)      # For data manipulation
library(tidyr)      # For data tidying
library(ggridges)   # For creating ridge plots
library(patchwork)  # For combining multiple plots into a single figure

# --- 2. Define File Paths ---
# Hardcoded paths are used as per user request.
# For better reproducibility, consider using relative paths with R Projects.
data_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
save_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"

# --- 3. Read and Prepare Data ---
effects_data <- read.csv(file.path(data_dir, "effect_sizes.csv"))

# --- 4. Define Plotting Aesthetics and Constants ---

# Define a consistent theme for all plots to ensure a uniform look
# with increased font sizes for better readability in publications.
publication_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = rel(1.4), face = "bold", hjust = 0.5),
    axis.title = element_text(size = rel(1.2), face = "bold"),
    axis.text = element_text(size = rel(1.0), color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Default angle
    legend.title = element_text(size = rel(1.2), face = "bold"),
    legend.text = element_text(size = rel(1.0)),
    legend.position = "top"
  )

# Define effect size thresholds for Cliff's Delta for visual reference on the plot.
# These typically represent small, medium, and large effects.
cliff_delta_thresholds <- c(-0.474, -0.33, -0.147, 0, 0.147, 0.33, 0.474)

# Define colors for genotypes to ensure consistency across plots.
genotype_colors <- c("G1" = "#2ecc71", "G2" = "#33d6d3")


# --- 5. Create Panel 1: Effect Size Distribution ---

p1_effect_distribution <- ggplot(effects_data, aes(x = effect_size, y = Tissue.type, fill = Genotype)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  geom_vline(xintercept = cliff_delta_thresholds,
             linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = genotype_colors) +
  scale_y_discrete(labels = c("L" = "Leaf", "R" = "Root")) +
  publication_theme +
  labs(
    title = "Effect Size Distribution by Tissue and Genotype",
    x = "Effect Size (Cliff's Delta, d)",
    y = "Tissue Type"
  )

# --- 6. Create Panel 2: Leaf-to-Root Response Ratio ---

# NOTE FOR PUBLICATION:
# The following data frame `ratio_data` is hardcoded. For full reproducibility,
# this data should be derived programmatically from the `effects_data` file.
# The logic would involve:
# 1. Defining effect size categories (e.g., |d| >= 0.474).
# 2. Counting the number of features within each category for each Genotype and Tissue.
# 3. Calculating the Leaf:Root count ratio for each Genotype and category.
#
# Example derivation:
# ratio_calculation <- effects_data %>%
#   mutate(abs_effect = abs(effect_size)) %>%
#   mutate(Category = case_when(
#     abs_effect >= 0.474 ~ ">= 0.474",
#     abs_effect >= 0.33  ~ "0.33-0.474",
#     TRUE ~ NA_character_
#   )) %>%
#   filter(!is.na(Category)) %>%
#   count(Genotype, Tissue.type, Category) %>%
#   tidyr::pivot_wider(names_from = Tissue.type, values_from = n, values_fill = 0) %>%
#   mutate(Ratio = ifelse(R == 0, 0, L / R)) # Avoid division by zero

ratio_data <- data.frame(
  Genotype = c("G1", "G1", "G2", "G2"),
  Category = factor(c("≥ 0.474", "0.33-0.474", "≥ 0.474", "0.33-0.474"),
                    levels = c("≥ 0.474", "0.33-0.474")),
  Ratio = c(10.2, 5.1, 0.4, 1.2)
)

p2_tissue_ratio <- ggplot(ratio_data, aes(x = Category, y = Ratio, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = genotype_colors) +
  publication_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate labels for better fit
  ) +
  labs(
    title = "Tissue Response Ratio",
    x = "Effect Size Category (|d|)",
    y = "Leaf : Root Ratio"
  )

# --- 7. Combine and Save the Final Plot ---

# Combine the two plots into a single figure using patchwork.
# The effect size distribution plot is given more width, and legends are collected.
combined_plot <- p1_effect_distribution + p2_tissue_ratio +
  plot_layout(ncol = 2, widths = c(3, 2), guides = 'collect') &
  theme(legend.position = 'top')

# Save the combined plot in both PDF (vector) and PNG (raster) formats.
ggsave(
  filename = file.path(save_dir, "effect_size_and_ratio.pdf"),
  plot = combined_plot,
  width = 12,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(save_dir, "effect_size_and_ratio.png"),
  plot = combined_plot,
  width = 12,
  height = 6,
  dpi = 300
)
