# -----------------------------------------------------------------------------
# Script: tissue_temporal.R
#
# Description:
# This script generates Figure 2, which visualizes the temporal dynamics and
# cross-tissue coordination of molecular features in response to osmotic stress.
# The figure consists of four panels:
#   A. Metabolic Trajectories: Line plot showing mean metabolite values over time.
#   B. Temporal Coordination: Line plot of cross-tissue correlations over time.
#   C. Tissue Coordination: Boxplot comparing metabolite values between tissues.
#   D. Response Metrics: Bar plot of mean molecular response by genotype and tissue.
#
# The final composite figure is saved in both PDF and TIFF formats.
#
# Author: [Your Name]
# Date: [Date]
# -----------------------------------------------------------------------------

# 1. SETUP
# -----------------------------------------------------------------------------
# Load required libraries
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)

# Set paths (as per user request, these are hardcoded)
out_dir <- "C:/Users/ms/Desktop/r/chem_data/metabo2/result/section3"
raw_data_path <- "C:/Users/ms/Desktop/r/chem_data/metabo2/Merged_VIP_hub_r_Path2.csv"

# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 2. DATA LOADING
# -----------------------------------------------------------------------------
# Read metabolite data
metabolite_data <- read.csv(raw_data_path)

# Prepare data: Add a 'Tissue' column with descriptive names
metabolite_data$Tissue <- ifelse(metabolite_data$Tissue.type == "L", "Leaf", "Root")

# 3. THEME AND PALETTE DEFINITION
# -----------------------------------------------------------------------------
# Define a color palette for consistency
genotype_colors <- c("G1" = "#2ecc71", "G2" = "#33d6d3")
tissue_colors <- c("Leaf" = "#27ae60", "Root" = "#33d6d3")

# Define a professional theme for plots (based on Nature's style)
nature_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# 4. PLOT GENERATION
# -----------------------------------------------------------------------------
# Panel A: Metabolic Trajectories
panel_a_data <- metabolite_data %>%
  group_by(Genotype, Tissue, Day) %>%
  summarise(
    mean_value = mean(Metabolite_Value, na.rm = TRUE),
    se_value = sd(Metabolite_Value, na.rm = TRUE) / sqrt(n()),
    n = n(),
    ci95 = qt(0.975, df = n - 1) * se_value,
    .groups = 'drop'
  )

metabolic_trajectories_plot <- ggplot(panel_a_data, aes(x = Day, y = mean_value, color = Genotype, group = interaction(Genotype, Tissue))) +
  geom_ribbon(aes(ymin = mean_value - ci95, ymax = mean_value + ci95, fill = Genotype), alpha = 0.2) +
  geom_line(size = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  facet_wrap(~Tissue, scales = "free_y") +
  labs(title = "Molecular Feature Trajectories",
       x = "Time (Days)",
       y = "Molecular Response (a.u.)") +
  scale_x_continuous(breaks = 1:3) +
  nature_theme

# Panel B: Temporal Coordination
# Note: Data for this panel is based on pre-calculated correlation values.
temporal_corr_data <- data.frame(
  Day = c(1, 2, 3, 1, 2, 3),
  Genotype = rep(c("G1", "G2"), each = 3),
  Correlation = c(0.546, 0.448, 0.350, 0.236, 0.262, 0.288),
  CI_Width = c(0.05, 0.05, 0.05, 0.04, 0.04, 0.04)
)

temporal_coordination_plot <- ggplot(temporal_corr_data, aes(x = Day, y = Correlation, color = Genotype, group = Genotype)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Correlation - CI_Width, ymax = Correlation + CI_Width), width = 0.1, alpha = 0.5) +
  scale_color_manual(values = genotype_colors) +
  scale_x_continuous(breaks = 1:3) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
  labs(title = "Temporal Coordination",
       x = "Time (Days)",
       y = "Cross-tissue Correlation (ρ)") +
  nature_theme

# Panel C: Tissue Coordination
tissue_coordination_plot <- ggplot(metabolite_data, aes(x = Tissue, y = Metabolite_Value, fill = Genotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), aes(color = Genotype), size = 1, alpha = 0.4) +
  scale_fill_manual(values = genotype_colors) +
  scale_color_manual(values = genotype_colors) +
  labs(title = "Tissue Coordination",
       y = "Molecular Response",
       x = "Tissue") +
  nature_theme

# Panel D: Response Metrics
panel_d_data <- metabolite_data %>%
  group_by(Tissue, Genotype) %>%
  summarise(
    mean_response = mean(abs(Metabolite_Value), na.rm = TRUE),
    sd_response = sd(abs(Metabolite_Value), na.rm = TRUE),
    n = n(),
    ci95 = qt(0.975, df = n - 1) * (sd_response / sqrt(n)),
    .groups = 'drop'
  )

response_metrics_plot <- ggplot(panel_d_data, aes(x = Genotype, y = mean_response, fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_response - ci95, ymax = mean_response + ci95), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = tissue_colors) +
  labs(title = "Response Metrics",
       y = "Mean Molecular Response",
       x = "Genotype") +
  nature_theme

# 5. PLOT ASSEMBLY
# -----------------------------------------------------------------------------
# Combine the four panels into a grid
combined_plot <- plot_grid(
  metabolic_trajectories_plot, temporal_coordination_plot,
  tissue_coordination_plot, response_metrics_plot,
  ncol = 2,
  align = 'hv',
  axis = 'tblr',
  rel_heights = c(1, 1)
)

# Add a main title to the combined plot
main_title <- ggdraw() + 
  draw_label(
    "Figure 2 | Temporal Dynamics & Cross-tissue Coordination in Response to Osmotic Stress",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 14
  ) +
  theme(plot.margin = margin(0, 0, 5, 7))

# Create the final plot with the title
final_plot <- plot_grid(
  main_title, combined_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

# 6. SAVE OUTPUT
# -----------------------------------------------------------------------------
# Save the final plot in multiple formats for publication
ggsave(file.path(out_dir, "Fig2_temporal_dynamics.pdf"),
       final_plot,
       width = 190,
       height = 160,
       units = "mm",
       dpi = 300)

ggsave(file.path(out_dir, "Fig2_temporal_dynamics.tiff"),
       final_plot,
       width = 190,
       height = 160,
       units = "mm",
       dpi = 300,
       compression = "lzw")

print("Plots saved successfully in PDF and TIFF formats.")
