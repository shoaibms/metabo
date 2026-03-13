# =============================================================================
# Figure 3 REDESIGN: 5-Panel Asymmetric Layout
# =============================================================================
#
# Layout:
#   Row 1: [A lollipop 1/3] [B violin+box 2/3]
#   Row 2: [C butterfly HERO -- full width]
#   Row 3: [D dumbbell 1/2] [E density overlay 1/2]
#
# Every panel uses a DIFFERENT visual grammar:
#   A - Lollipop:    Hub persistence (Jaccard similarity)
#   B - Violin+Box:  VIP temporal slope distributions
#   C - Butterfly:   Dynamism (left) vs Resilience (right) -- HERO
#   D - Dumbbell:    Mean |RMI| tissue contrast per genotype
#   E - Density:     Observed vs null effect-size distributions
#
# Colour scheme: #3CB371 (G1) / #2AA9B0 (G2) matching fig_2_a_b_c.R
#
# Author: Plant Systems Biology Lab
# Date: 2025 (Revision)
# =============================================================================

# 1. ENVIRONMENT SETUP
# =============================================================================
rm(list = ls())
gc()

required_packages <- c(
  "tidyverse", "ggplot2", "cowplot", "viridis",
  "RColorBrewer", "scales", "grid", "gridExtra",
  "effsize", "boot", "reshape2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 2. FILE PATHS
# >> USER CONFIGURATION — Update all paths below for your local environment. <<
# =============================================================================
data_paths <- list(
  leaf_data = "C:/Users/USER/Desktop/data_chem_3_10/data/data/n_p_l.csv",
  root_data = "C:/Users/USER/Desktop/data_chem_3_10/data/data/n_p_r.csv",
  vip_data  = "C:/Users/USER/Desktop/data_chem_3_10/output/results/spearman/VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv",
  temporal  = "C:/Users/USER/Desktop/data_chem_3_10/output/results/initial_stat/initial_stat6e/temporal_analysis.csv"
)

hub_paths <- list(
  leaf_g1 = "C:/Users/USER/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/leaf_g1_hub_metabolites.csv",
  leaf_g2 = "C:/Users/USER/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/leaf_g2_hub_metabolites.csv",
  root_g1 = "C:/Users/USER/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/root_g1_hub_metabolites.csv",
  root_g2 = "C:/Users/USER/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/root_g2_hub_metabolites.csv"
)

output_dir <- "C:/Users/USER/Desktop/data_chem_3_10/output/revision_pack/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. DESIGN SYSTEM
# =============================================================================
COL_G1       <- "#3CB371"    # Medium sea green (tolerant)
COL_G2       <- "#2AA9B0"    # Teal (susceptible)
COL_G1_DARK  <- "#2E8B57"    # Sea green (darker G1)
COL_G2_DARK  <- "#1B6F75"    # Dark teal (darker G2)
COL_ACCENT   <- "#DC6B4A"    # Warm coral (annotations, reference lines)
COL_NULL     <- "#F5E6A3"    # Warm light yellow (null distributions -- vibrancy)
COL_BG       <- "white"
COL_GRID     <- "grey92"
COL_BORDER   <- "black"
COL_ANNOT    <- "#34495e"    # Dark slate (annotation text)

palette_geno   <- c("G1" = COL_G1, "G2" = COL_G2)
palette_tissue <- c("Leaf" = COL_G1, "Root" = COL_G2)

# Shared theme
base_theme <- function(base_sz = 14) {
  theme_minimal(base_size = base_sz) +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      axis.title       = element_text(size = 16, face = "plain"),
      axis.title.y     = element_text(size = 16, margin = margin(r = 8)),
      axis.text        = element_text(size = 16, color = "black"),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.position  = "top",
      legend.justification = "left",
      legend.box.spacing   = unit(0.2, "cm"),
      legend.key.size      = unit(0.5, "cm"),
      legend.margin        = margin(2, 0, 2, 0),
      panel.grid.major = element_line(color = COL_GRID, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = COL_BORDER, fill = NA, linewidth = 0.7),
      strip.text       = element_text(size = 16, face = "bold"),
      strip.background = element_rect(fill = "grey96", color = COL_BORDER, linewidth = 0.5),
      plot.margin      = margin(6, 3, 2, 3),
      plot.tag          = element_text(size = 32, face = "plain"),
      plot.tag.position = "topleft"
    )
}

# 4. DATA LOADING
# =============================================================================
cat("Loading data...\n")

safe_read <- function(path, desc) {
  if (!file.exists(path)) { cat("Warning: Missing:", path, "\n"); return(NULL) }
  tryCatch({
    d <- read.csv(path, stringsAsFactors = FALSE)
    cat("[OK]", desc, ":", nrow(d), "x", ncol(d), "\n")
    d
  }, error = function(e) { cat("[FAIL]", desc, "\n"); NULL })
}

leaf_data <- safe_read(data_paths$leaf_data, "leaf")
root_data <- safe_read(data_paths$root_data, "root")
vip_data  <- safe_read(data_paths$vip_data, "VIP")
temp_data <- safe_read(data_paths$temporal, "temporal")

hub_data <- list()
for (nm in names(hub_paths)) hub_data[[nm]] <- safe_read(hub_paths[[nm]], nm)
hub_data <- hub_data[!sapply(hub_data, is.null)]

get_metab_cols <- function(d) {
  if (is.null(d)) return(character(0))
  grep("^(N_Cluster_|P_Cluster_)", names(d), value = TRUE)
}

shared_metabolites <- intersect(get_metab_cols(leaf_data), get_metab_cols(root_data))
cat("[OK]", length(shared_metabolites), "shared metabolites\n")

# 5. ANALYSIS FUNCTIONS
# =============================================================================

# --- Hub persistence ---
calc_hub_persistence <- function(hub_data) {
  cat("Hub persistence...\n")
  if (length(hub_data) == 0) return(NULL)
  
  hub_sets <- list()
  for (cond in names(hub_data)) {
    hd <- hub_data[[cond]]
    if (!is.null(hd) && nrow(hd) > 0) {
      hub_sets[[cond]] <- if ("Metabolite" %in% names(hd)) hd$Metabolite else hd[, 1]
    }
  }
  
  rows <- list()
  # Genotype comparisons within tissues
  for (tis in c("leaf", "root")) {
    g1k <- paste0(tis, "_g1"); g2k <- paste0(tis, "_g2")
    if (g1k %in% names(hub_sets) && g2k %in% names(hub_sets)) {
      g1h <- hub_sets[[g1k]]; g2h <- hub_sets[[g2k]]
      sh <- intersect(g1h, g2h)
      rows[[length(rows) + 1]] <- data.frame(
        Label = paste(tools::toTitleCase(tis), "G1 vs G2"),
        Type = "Genotype", Jaccard = length(sh) / length(union(g1h, g2h)),
        stringsAsFactors = FALSE)
    }
  }
  # Tissue comparisons within genotypes
  for (gen in c("g1", "g2")) {
    lk <- paste0("leaf_", gen); rk <- paste0("root_", gen)
    if (lk %in% names(hub_sets) && rk %in% names(hub_sets)) {
      lh <- hub_sets[[lk]]; rh <- hub_sets[[rk]]
      sh <- intersect(lh, rh)
      rows[[length(rows) + 1]] <- data.frame(
        Label = paste(toupper(gen), "Leaf vs Root"),
        Type = "Tissue", Jaccard = length(sh) / length(union(lh, rh)),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# --- VIP tracking ---
calc_vip_tracking <- function(leaf_data, root_data, vip_data, shared_metabolites) {
  cat("VIP tracking...\n")
  if (is.null(vip_data) || nrow(vip_data) == 0) return(NULL)
  vip_in <- intersect(vip_data$Metabolite, shared_metabolites)
  if (length(vip_in) == 0) return(NULL)
  
  common <- intersect(names(leaf_data), names(root_data))
  req_cols <- c("Day", "Genotype", "Treatment")
  met_cols <- intersect(shared_metabolites, common)
  all_cols <- unique(c(req_cols, met_cols))
  
  comb <- rbind(
    leaf_data %>% select(all_of(all_cols)) %>% mutate(Tissue = "Leaf"),
    root_data %>% select(all_of(all_cols)) %>% mutate(Tissue = "Root")
  )
  
  res <- list()
  for (m in vip_in) {
    for (g in c("G1", "G2")) {
      for (tis in c("Leaf", "Root")) {
        for (tr in c(0, 1)) {
          sub <- comb %>% filter(Genotype == g, Tissue == tis, Treatment == tr) %>%
            select(Day, all_of(m))
          if (nrow(sub) == 0 || !(m %in% names(sub))) next
          
          ts <- sub %>% group_by(Day) %>%
            summarise(mv = mean(.data[[m]], na.rm = TRUE),
                      n = sum(!is.na(.data[[m]])), .groups = "drop") %>%
            filter(n > 0)
          if (nrow(ts) < 2) next
          
          slope <- NA
          if (nrow(ts) >= 3) {
            tryCatch({
              lmr <- lm(mv ~ Day, data = ts); slope <- coef(lmr)[2]
            }, error = function(e) { slope <<- NA })
          } else {
            slope <- (ts$mv[2] - ts$mv[1]) / (ts$Day[2] - ts$Day[1])
          }
          
          key <- paste(m, g, tis, tr, sep = "_")
          res[[key]] <- data.frame(
            Metabolite = m, Genotype = g, Tissue = tis, Treatment = tr,
            Temporal_slope = ifelse(is.na(slope), 0, slope),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  do.call(rbind, res)
}

# --- Feature changes (dynamism) ---
calc_feature_changes <- function(leaf_data, root_data, shared_metabolites) {
  cat("Feature changes...\n")
  
  test_trends <- function(data, tissue_name) {
    rows <- list()
    for (g in c("G1", "G2")) {
      sub <- data %>% filter(Genotype == g, Treatment == 1) %>%
        select(Day, all_of(shared_metabolites))
      if (nrow(sub) == 0) next
      sig_count <- sum(sapply(shared_metabolites, function(m) {
        md <- sub %>% select(Day, all_of(m)) %>% filter(!is.na(.data[[m]]))
        if (nrow(md) >= 3) {
          lmr <- tryCatch(lm(get(m) ~ Day, data = md), error = function(e) NULL)
          if (!is.null(lmr)) return(summary(lmr)$coefficients[2, 4] < 0.05)
        }
        FALSE
      }), na.rm = TRUE)
      rows[[g]] <- data.frame(Genotype = g, Tissue = tissue_name,
                              Dynamic_features = sig_count, stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  }
  
  rbind(test_trends(leaf_data, "Leaf"), test_trends(root_data, "Root"))
}

# --- Feature stability (resilience) ---
calc_feature_stability <- function(leaf_data, root_data, shared_metabolites) {
  cat("Feature stability...\n")
  
  calc_res <- function(data, tissue_name) {
    rows <- list()
    for (g in c("G1", "G2")) {
      ctrl <- data %>% filter(Genotype == g, Treatment == 0) %>%
        select(all_of(shared_metabolites)) %>% summarise_all(mean, na.rm = TRUE)
      strs <- data %>% filter(Genotype == g, Treatment == 1,
                              Day == max(Day, na.rm = TRUE)) %>%
        select(all_of(shared_metabolites)) %>% summarise_all(mean, na.rm = TRUE)
      if (nrow(ctrl) == 0 || nrow(strs) == 0) next
      rc <- 0; tc <- 0
      for (m in shared_metabolites) {
        cv <- ctrl[[m]]; sv <- strs[[m]]
        if (!is.na(cv) && !is.na(sv) && cv != 0) {
          if (abs(sv - cv) / abs(cv) < 0.2) rc <- rc + 1
          tc <- tc + 1
        }
      }
      rows[[g]] <- data.frame(Genotype = g, Tissue = tissue_name,
                              Resilient_pct = rc / tc * 100, stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  }
  
  rbind(calc_res(leaf_data, "Leaf"), calc_res(root_data, "Root"))
}

# --- Mean response ---
calc_mean_responses <- function(leaf_data, root_data, shared_metabolites) {
  cat("Mean responses...\n")
  
  proc <- function(data, tis) {
    data %>% filter(Treatment == 1) %>%
      select(Genotype, all_of(shared_metabolites)) %>%
      pivot_longer(cols = all_of(shared_metabolites), names_to = "M", values_to = "V") %>%
      group_by(Genotype) %>%
      summarise(Mean_RMI = mean(abs(V), na.rm = TRUE),
                SD = sd(abs(V), na.rm = TRUE),
                N = sum(!is.na(V)), .groups = "drop") %>%
      mutate(Tissue = tis, SE = SD / sqrt(N),
             CI_lo = Mean_RMI - 1.96 * SE, CI_hi = Mean_RMI + 1.96 * SE)
  }
  
  bind_rows(proc(leaf_data, "Leaf"), proc(root_data, "Root"))
}

# --- Null distributions ---
calc_null_dist <- function(leaf_data, root_data, shared_metabolites, n_perm = 500) {
  cat("Null distributions (", n_perm, "permutations)...\n")
  samp_met <- sample(shared_metabolites, min(50, length(shared_metabolites)))
  null_res <- list()
  
  for (i in 1:n_perm) {
    if (i %% 100 == 0) cat("  ", i, "/", n_perm, "\n")
    for (g in c("G1", "G2")) {
      for (td in list(list(d = leaf_data, t = "Leaf"), list(d = root_data, t = "Root"))) {
        gd <- td$d %>% filter(Genotype == g) %>% select(Treatment, all_of(samp_met))
        if (nrow(gd) <= 10) next
        pd <- gd; pd$Treatment <- sample(gd$Treatment)
        m <- samp_met[1]
        c_v <- pd[pd$Treatment == 0, m]; t_v <- pd[pd$Treatment == 1, m]
        cc <- c_v[!is.na(c_v)]; tc <- t_v[!is.na(t_v)]
        if (length(cc) >= 3 && length(tc) >= 3) {
          tryCatch({
            cr <- suppressWarnings(cliff.delta(tc, cc))
            null_res[[paste(i, g, td$t, sep = "_")]] <- data.frame(
              Genotype = g, Tissue = td$t, Null_ES = cr$estimate)
          }, error = function(e) {})
        }
      }
    }
  }
  do.call(rbind, null_res)
}

# --- Effect sizes (observed) ---
calc_effect_sizes <- function(leaf_data, root_data, shared_metabolites) {
  cat("Observed effect sizes...\n")
  res <- list()
  for (g in c("G1", "G2")) {
    for (td in list(list(d = leaf_data, t = "Leaf"), list(d = root_data, t = "Root"))) {
      ctrl <- td$d %>% filter(Genotype == g, Treatment == 0) %>% select(all_of(shared_metabolites))
      trt  <- td$d %>% filter(Genotype == g, Treatment == 1) %>% select(all_of(shared_metabolites))
      if (nrow(ctrl) == 0 || nrow(trt) == 0) next
      for (m in shared_metabolites) {
        cc <- ctrl[[m]][!is.na(ctrl[[m]])]; tc <- trt[[m]][!is.na(trt[[m]])]
        if (length(cc) >= 3 && length(tc) >= 3) {
          tryCatch({
            cr <- suppressWarnings(cliff.delta(tc, cc))
            res[[paste(g, td$t, m, sep = "_")]] <- data.frame(
              Genotype = g, Tissue = td$t, Metabolite = m,
              Effect_Size = cr$estimate)
          }, error = function(e) {})
        }
      }
    }
  }
  do.call(rbind, res)
}

# 6. RUN ALL ANALYSES
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("RUNNING ANALYSES\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

hub_data_res     <- calc_hub_persistence(hub_data)
vip_tracking     <- calc_vip_tracking(leaf_data, root_data, vip_data, shared_metabolites)
feature_dynamics <- calc_feature_changes(leaf_data, root_data, shared_metabolites)
feature_stab     <- calc_feature_stability(leaf_data, root_data, shared_metabolites)
mean_resp        <- calc_mean_responses(leaf_data, root_data, shared_metabolites)
null_dist        <- calc_null_dist(leaf_data, root_data, shared_metabolites, 500)
obs_effects      <- calc_effect_sizes(leaf_data, root_data, shared_metabolites)

cat("[OK] All analyses complete\n\n")

# 7. PANEL CREATION
# =============================================================================

# ---- Panel A: LOLLIPOP (Hub Persistence) ----
create_panel_a <- function(hub_df) {
  cat("Panel A: Lollipop...\n")
  if (is.null(hub_df)) return(ggplot() + theme_void() + labs(tag = "A"))
  
  hub_df <- hub_df %>% arrange(Jaccard)
  
  # Ordered factor for y-axis
  hub_df$Label <- factor(hub_df$Label, levels = hub_df$Label)
  
  ggplot(hub_df, aes(x = Jaccard, y = Label)) +
    # Stem
    geom_segment(aes(x = 0, xend = Jaccard, y = Label, yend = Label, color = Type),
                 linewidth = 1.2, lineend = "round") +
    # Dot
    geom_point(aes(color = Type), size = 5) +
    # Value label
    geom_text(aes(label = sprintf("%.3f", Jaccard)),
              hjust = -0.4, fontface = "bold", size = 5, color = COL_ANNOT) +
    scale_color_manual(
      values = c("Genotype" = COL_G1, "Tissue" = COL_G2),
      name = NULL,
      labels = c("Genotype" = "Between genotypes", "Tissue" = "Between tissues")
    ) +
    scale_x_continuous(limits = c(0, max(hub_df$Jaccard) * 1.55),
                       labels = number_format(accuracy = 0.1)) +
    labs(x = "Jaccard Index", y = NULL, tag = "A") +
    base_theme() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 16)
    )
}

# ---- Panel B: VIOLIN+BOX (VIP Temporal Dynamics) ----
create_panel_b <- function(vip_df) {
  cat("Panel B: Violin+Box...\n")
  if (is.null(vip_df)) return(ggplot() + theme_void() + labs(tag = "B"))
  
  stress_vip <- vip_df %>% filter(Treatment == 1, !is.na(Temporal_slope))
  if (nrow(stress_vip) == 0) return(ggplot() + theme_void() + labs(tag = "B"))
  
  # IQR-based zoom
  q1 <- quantile(stress_vip$Temporal_slope, 0.25, na.rm = TRUE)
  q3 <- quantile(stress_vip$Temporal_slope, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  y_lo <- q1 - 2.5 * iqr; y_hi <- q3 + 2.5 * iqr
  
  ggplot(stress_vip, aes(x = Tissue, y = Temporal_slope,
                         fill = Genotype, color = Genotype)) +
    geom_violin(alpha = 0.45, position = position_dodge(0.85),
                linewidth = 0.5, scale = "width") +
    geom_boxplot(width = 0.18, position = position_dodge(0.85),
                 outlier.shape = NA, alpha = 0.9, linewidth = 0.5,
                 color = "grey30") +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = COL_ACCENT, alpha = 0.7, linewidth = 0.6) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    scale_fill_manual(values = palette_geno, name = "Genotype",
                      labels = c("G1" = "G1 (tolerant)", "G2" = "G2 (susceptible)")) +
    scale_color_manual(values = palette_geno, name = "Genotype",
                       labels = c("G1" = "G1 (tolerant)", "G2" = "G2 (susceptible)")) +
    labs(x = "Tissue", y = "Temporal Slope", tag = "B") +
    base_theme() +
    theme(legend.position = "top")
}

# ---- Panel C: BUTTERFLY / DIVERGING BAR (Hero -- Dynamism vs Resilience) ----
create_panel_c <- function(dynamics_df, stability_df) {
  cat("Panel C: Butterfly (Hero)...\n")
  if (is.null(dynamics_df) || is.null(stability_df)) {
    return(ggplot() + theme_void() + labs(tag = "C"))
  }
  
  # Merge dynamism and resilience data
  combo <- dynamics_df %>%
    inner_join(stability_df, by = c("Genotype", "Tissue")) %>%
    mutate(
      Row_Label = paste0(Genotype, " ", Tissue),
      # Negate dynamism for left side
      Dynamic_neg = -Dynamic_features
    )
  
  # Order: G1 Leaf, G1 Root, G2 Leaf, G2 Root (top to bottom)
  row_order <- c("G1 Leaf", "G1 Root", "G2 Leaf", "G2 Root")
  combo$Row_Label <- factor(combo$Row_Label, levels = rev(row_order))
  
  # Assign colours per genotype
  combo$Fill_Col <- ifelse(combo$Genotype == "G1", COL_G1, COL_G2)
  
  # Compute deltas for annotations
  g1_dyn_delta <- combo$Dynamic_features[combo$Row_Label == "G1 Root"] -
    combo$Dynamic_features[combo$Row_Label == "G1 Leaf"]
  g1_res_delta <- combo$Resilient_pct[combo$Row_Label == "G1 Root"] -
    combo$Resilient_pct[combo$Row_Label == "G1 Leaf"]
  
  # Max values for scaling
  max_dyn <- max(combo$Dynamic_features) * 1.25
  max_res <- max(combo$Resilient_pct) * 1.35
  
  # --- Left panel: Dynamism (features with sig trends) ---
  p_left <- ggplot(combo, aes(x = Dynamic_neg, y = Row_Label)) +
    geom_col(aes(fill = Genotype), width = 0.65, alpha = 0.85) +
    geom_text(aes(label = Dynamic_features, x = Dynamic_neg - max_dyn * 0.04),
              hjust = 1, fontface = "bold", size = 5.5, color = COL_ANNOT) +
    scale_fill_manual(values = palette_geno, guide = "none") +
    scale_x_continuous(
      limits = c(-max_dyn, 0),
      labels = function(x) abs(x),
      expand = expansion(mult = c(0.05, 0))
    ) +
    labs(x = "Dynamic features (sig. temporal trends)", y = NULL,
         title = "Temporal Dynamism") +
    base_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15, margin = margin(b = 10)),
      plot.margin = margin(2, 0, 2, 6)
    )
  
  # --- Centre strip: Row labels ---
  p_centre <- ggplot(combo, aes(x = 0, y = Row_Label)) +
    geom_text(aes(label = Row_Label, color = Genotype),
              size = 5.5, fontface = "bold") +
    scale_color_manual(values = palette_geno, guide = "none") +
    scale_x_continuous(limits = c(-1, 1)) +
    theme_void() +
    theme(plot.margin = margin(2, 0, 2, 0))
  
  # --- Right panel: Resilience (% stable features) ---
  p_right <- ggplot(combo, aes(x = Resilient_pct, y = Row_Label)) +
    geom_col(aes(fill = Genotype), width = 0.65, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", Resilient_pct),
                  x = Resilient_pct + max_res * 0.03),
              hjust = 0, fontface = "bold", size = 5.5, color = COL_ANNOT) +
    scale_fill_manual(values = palette_geno, guide = "none") +
    scale_x_continuous(
      limits = c(0, max_res),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = "Resilient features (%)", y = NULL,
         title = "Metabolic Resilience") +
    base_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15, margin = margin(b = 10)),
      plot.margin = margin(2, 6, 2, 0)
    )
  
  # --- Assemble butterfly ---
  butterfly_body <- plot_grid(
    p_left, p_centre, p_right,
    ncol = 3, align = "h", axis = "tb",
    rel_widths = c(1, 0.35, 1)
  )
  
  # Tag
  tag_grob <- ggdraw() +
    draw_label("C", fontface = "plain", size = 32, hjust = 0, x = 0.005, y = 0.7)
  
  # Combine
  plot_grid(
    tag_grob,
    butterfly_body,
    ncol = 1, rel_heights = c(0.07, 1)
  )
}

# ---- Panel D: DUMBBELL (Mean |RMI| tissue contrast) ----
create_panel_d <- function(mean_resp) {
  cat("Panel D: Dumbbell...\n")
  if (is.null(mean_resp)) return(ggplot() + theme_void() + labs(tag = "D"))
  
  # Pivot wider for dumbbell
  wide <- mean_resp %>%
    select(Genotype, Tissue, Mean_RMI, CI_lo, CI_hi) %>%
    pivot_wider(names_from = Tissue,
                values_from = c(Mean_RMI, CI_lo, CI_hi))
  
  wide$Genotype_label <- factor(
    ifelse(wide$Genotype == "G1", "G1 (tolerant)", "G2 (susceptible)"),
    levels = c("G2 (susceptible)", "G1 (tolerant)")
  )
  
  ggplot(wide) +
    # Connecting segment (the "dumbbell bar")
    geom_segment(aes(x = Mean_RMI_Leaf, xend = Mean_RMI_Root,
                     y = Genotype_label, yend = Genotype_label),
                 linewidth = 1.8, color = "grey70", lineend = "round") +
    # Error bars
    geom_errorbarh(aes(xmin = CI_lo_Leaf, xmax = CI_hi_Leaf,
                       y = Genotype_label),
                   height = 0.15, linewidth = 0.5, color = COL_G1_DARK) +
    geom_errorbarh(aes(xmin = CI_lo_Root, xmax = CI_hi_Root,
                       y = Genotype_label),
                   height = 0.15, linewidth = 0.5, color = COL_G2_DARK) +
    # Dots
    geom_point(aes(x = Mean_RMI_Leaf, y = Genotype_label),
               size = 6, color = COL_G1, shape = 16) +
    geom_point(aes(x = Mean_RMI_Root, y = Genotype_label),
               size = 6, color = COL_G2, shape = 17) +
    # Value labels
    geom_text(aes(x = Mean_RMI_Leaf, y = Genotype_label,
                  label = round(Mean_RMI_Leaf, 0)),
              vjust = -1.3, size = 5, fontface = "bold", color = COL_G1_DARK) +
    geom_text(aes(x = Mean_RMI_Root, y = Genotype_label,
                  label = round(Mean_RMI_Root, 0)),
              vjust = -1.3, size = 5, fontface = "bold", color = COL_G2_DARK) +
    # Tissue delta annotation
    geom_text(aes(x = (Mean_RMI_Leaf + Mean_RMI_Root) / 2,
                  y = Genotype_label,
                  label = paste0("D = ", round(abs(Mean_RMI_Root - Mean_RMI_Leaf), 0))),
              vjust = 2.5, size = 4.5, fontface = "italic", color = COL_ANNOT) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    labs(x = "Mean |RMI| (a.u.)", y = NULL, tag = "D") +
    base_theme() +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "top"
    ) +
    # Manual legend for shapes
    annotate("point", x = -Inf, y = -Inf, shape = 16, color = COL_G1, size = 0) +
    annotate("point", x = -Inf, y = -Inf, shape = 17, color = COL_G2, size = 0)
}

# Add a manual legend to Panel D using a small helper
add_dumbbell_legend <- function(p) {
  p + labs(caption = NULL) +
    # Use annotation_custom or just rely on the subtitle being self-explanatory
    # Add a legend via cowplot overlay
    theme(plot.subtitle = element_text(size = 17, color = "grey35"))
}

# ---- Panel E: DENSITY OVERLAY (Observed vs Null) ----
create_panel_e <- function(obs_df, null_df) {
  cat("Panel E: Density overlay...\n")
  if (is.null(obs_df) || is.null(null_df)) {
    return(ggplot() + theme_void() + labs(tag = "E"))
  }
  
  # Prepare observed
  obs_long <- obs_df %>%
    select(Genotype, Tissue, Effect_Size) %>%
    mutate(Source = "Observed")
  
  # Prepare null
  null_long <- null_df %>%
    select(Genotype, Tissue, Effect_Size = Null_ES) %>%
    mutate(Source = "Null (permuted)")
  
  # Combine
  all_data <- bind_rows(obs_long, null_long) %>%
    mutate(Facet_Label = paste(Tissue, Genotype, sep = " - "))
  
  ggplot(all_data, aes(x = abs(Effect_Size))) +
    # Null as grey fill
    geom_density(data = all_data %>% filter(Source == "Null (permuted)"),
                 fill = COL_NULL, color = "#D4C572", alpha = 0.6, linewidth = 0.4) +
    # Observed as coloured line + subtle fill
    geom_density(data = all_data %>% filter(Source == "Observed"),
                 aes(fill = Tissue, color = Tissue),
                 alpha = 0.35, linewidth = 0.8) +
    facet_wrap(~ Facet_Label, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = palette_tissue, name = "Observed") +
    scale_color_manual(values = c("Leaf" = COL_G1_DARK, "Root" = COL_G2_DARK),
                       guide = "none") +
    scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
    labs(x = "|Cliff's d|", y = "Density", tag = "E") +
    base_theme() +
    theme(
      strip.text = element_text(size = 16, face = "bold"),
      legend.position = "top"
    )
}

# 8. ASSEMBLE FIGURE
# =============================================================================
cat("Assembling panels...\n")

pa <- create_panel_a(hub_data_res)
pb <- create_panel_b(vip_tracking)
pc <- create_panel_c(feature_dynamics, feature_stab)

# Nudge panel B slightly to the right (~6-8 spacebar equivalent)
# Use left plot margin so panel height/aspect is preserved.
pb_shifted <- pb + theme(plot.margin = margin(6, 3, 2, 32))

# Add a small legend strip to Panel D for Leaf (circle) vs Root (triangle)
pd_raw <- create_panel_d(mean_resp)
# Create a tiny tissue legend
tissue_legend_df <- data.frame(
  x = c(1, 2), y = c(1, 1),
  Tissue = c("Leaf", "Root")
)
tissue_legend_plot <- ggplot(tissue_legend_df, aes(x, y, shape = Tissue, color = Tissue)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c("Leaf" = 16, "Root" = 17), name = "Tissue") +
  scale_color_manual(values = palette_tissue, name = "Tissue") +
  theme_void() +
  theme(legend.position = "top",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17, face = "bold"))
# Extract just the legend
tissue_legend <- get_legend(tissue_legend_plot)

# Add the legend above Panel D
pd <- plot_grid(tissue_legend, pd_raw, ncol = 1, rel_heights = c(0.08, 1))

pe <- create_panel_e(obs_effects, null_dist)

# --- ASYMMETRIC LAYOUT ---
cat("Building asymmetric layout...\n")

# Row 1: A (1/3) + B (2/3)
row1 <- plot_grid(pa, pb_shifted, ncol = 2, align = "h", axis = "tb",
                  rel_widths = c(1, 1.8))

# Row 2: C (full width hero)
row2 <- pc

# Row 3: D (1/2) + E (1/2)
row3 <- plot_grid(pd, pe, ncol = 2, align = "h", axis = "tb",
                  rel_widths = c(0.85, 1.15))

# Final assembly
final_figure <- plot_grid(
  row1, row2, row3,
  ncol = 1,
  rel_heights = c(1, 0.85, 1.3)
)

# 9. SAVE
# =============================================================================
cat("Saving...\n")

tryCatch({
  ggsave(
    filename = file.path(output_dir, "Fig3_Redesign_5Panel.png"),
    plot = final_figure,
    width = 16.5, height = 12.8,
    dpi = 300, units = "in", bg = "white"
  )
  cat("[OK] Saved: Fig3_Redesign_5Panel.png\n")
}, error = function(e) cat("[FAIL]", e$message, "\n"))

# Essential CSVs
if (!is.null(hub_data_res))
  write.csv(hub_data_res, file.path(output_dir, "fig3_redesign_hub_persistence.csv"), row.names = FALSE)
if (!is.null(feature_dynamics))
  write.csv(feature_dynamics, file.path(output_dir, "fig3_redesign_dynamics.csv"), row.names = FALSE)
if (!is.null(feature_stab))
  write.csv(feature_stab, file.path(output_dir, "fig3_redesign_stability.csv"), row.names = FALSE)

cat("\n[DONE] Figure 3 redesign complete\n")
cat("Output:", file.path(output_dir, "Fig3_Redesign_5Panel.png"), "\n")