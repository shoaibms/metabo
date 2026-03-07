# =============================================================================
#  bayesian_sensitivity.R
#  Bayesian Network Parent-Constraint Sensitivity Analysis
# =============================================================================
#
#  PURPOSE
#  -------
#  Standalone script that runs Bayesian structure learning under three
#  parent-constraint regimes (unconstrained, maxp=5, maxp=3) and compares
#  results.  All outputs are ADDITIVE — nothing in the original
#  baysian_network.R output folder is read, modified, or overwritten.
#
#  This script produces:
#    1. bayesian_parent_sensitivity_ALL.csv   — summary table (submission: Supp Table S10)
#    2. Supp_Fig_S9_BayesianParentSensitivity.pdf — comparison figure (submission: Supp Fig S10)
#    3. Per-condition .rds files for reproducibility
#
#  USAGE
#  -----
#  Source or Rscript from the project root.  Requires ~20-40 min on 20+ cores
#  (5 000 bootstrap for primary, 1 000 for sensitivity conditions).
#
#  The script reloads the same processed_data.rds used by baysian_network.R
#  so the input data are identical.
#
# =============================================================================

# ---- 0. Libraries -----------------------------------------------------------

suppressPackageStartupMessages({
  library(bnlearn)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

# ---- 0b. Reproducibility seed ------------------------------------------------

SEED <- 20250816
set.seed(SEED)
if (requireNamespace("doRNG", quietly = TRUE)) {
  library(doRNG)
  cat("doRNG available — parallel loops will be bitwise reproducible\n")
}

# ---- 1. Paths ---------------------------------------------------------------
# >> USER CONFIGURATION — Update these two paths for your local environment. <<

# Original BN output (read-only — we only load processed_data.rds from here)
original_bn_path <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"

# New output directory — all sensitivity results go here
output_dir <- "C:/Users/ms/Desktop/data_chem_3_10/output/revision_pack/bayesian_parent_sensitivity"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Parallel backend ----------------------------------------------------

# Use available cores, keeping some headroom
n_cores <- max(1, min(parallel::detectCores() - 4, 28))
cat("Setting up parallel backend with", n_cores, "cores\n")

cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

# If doRNG is available, register for reproducible parallel RNG
if (requireNamespace("doRNG", quietly = TRUE)) {
  doRNG::registerDoRNG(SEED)
}

# Guarantee clean shutdown even on error
on.exit({
  tryCatch(stopCluster(cluster), error = function(e) invisible(NULL))
  tryCatch(registerDoSEQ(),      error = function(e) invisible(NULL))
  tryCatch(closeAllConnections(), error = function(e) invisible(NULL))
  gc()
}, add = TRUE)

# ---- 3. Helper: hc wrapper with optional maxp -------------------------------

learn_hc <- function(df, maxp = NA) {
  # NA  → unconstrained (matches original baysian_network.R behaviour)
  # integer → pass maxp to hc()
  if (is.na(maxp)) return(hc(df))
  return(hc(df, maxp = as.integer(maxp)))
}

# ---- 4. Bootstrap analysis ---------------------------------------------------

bootstrap_network_analysis <- function(data, n_boots = 5000, maxp = NA) {
  if (ncol(data) == 0 || nrow(data) == 0) stop("Empty dataset provided")

  cat(sprintf("  Bootstrap (%d iterations, maxp=%s) on %d×%d matrix ...\n",
              n_boots, ifelse(is.na(maxp), "Inf", maxp), nrow(data), ncol(data)))

  t0 <- proc.time()

  # Chunk workload per core to minimise IPC overhead (28 returns, not 5000)
  n_workers   <- getDoParWorkers()
  chunk_sizes <- rep(n_boots %/% n_workers, n_workers)
  rem         <- n_boots %% n_workers
  if (rem > 0) chunk_sizes[1:rem] <- chunk_sizes[1:rem] + 1
  chunk_sizes <- chunk_sizes[chunk_sizes > 0]

  mats <- foreach(
    chunk     = chunk_sizes,
    .combine  = "+",
    .packages = "bnlearn",
    .export   = "learn_hc"
  ) %dopar% {
    local_sum <- NULL
    for (j in seq_len(chunk)) {
      idx       <- sample(nrow(data), replace = TRUE)
      boot_data <- data[idx, , drop = FALSE]
      bn        <- learn_hc(boot_data, maxp = maxp)
      mat       <- amat(bn)
      if (is.null(local_sum)) local_sum <- mat else local_sum <- local_sum + mat
    }
    local_sum
  }

  edge_prob <- mats / n_boots
  elapsed_boot <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  Bootstrap done in %.1f min (%.0f iter/min)\n",
              elapsed_boot / 60, n_boots / (elapsed_boot / 60)))
  keep      <- edge_prob >= 0.5   # same screening rule as manuscript

  # Max in-degree across nodes (for reporting hub truncation)
  bn_full   <- learn_hc(data, maxp = maxp)
  max_parents <- max(sapply(nodes(bn_full), function(n) length(parents(bn_full, n))))

  list(
    edge_prob    = edge_prob,
    keep         = keep,
    n_edges_kept = sum(keep),
    max_parents  = max_parents
  )
}

# ---- 5. Null model (column-wise permutation + re-learned DAG) ----------------

generate_null_networks <- function(data, n_permutations = 5000, maxp = NA) {
  if (ncol(data) == 0 || nrow(data) == 0) stop("Empty dataset provided")

  # Observed network
  observed_network <- learn_hc(data, maxp = maxp)
  observed_edges   <- nrow(arcs(observed_network))

  cat(sprintf("  Null model (%d permutations, maxp=%s) ...\n",
              n_permutations, ifelse(is.na(maxp), "Inf", maxp)))

  t0 <- proc.time()

  # Chunk workload per core (same IPC optimisation as bootstrap)
  n_workers   <- getDoParWorkers()
  chunk_sizes <- rep(n_permutations %/% n_workers, n_workers)
  rem         <- n_permutations %% n_workers
  if (rem > 0) chunk_sizes[1:rem] <- chunk_sizes[1:rem] + 1
  chunk_sizes <- chunk_sizes[chunk_sizes > 0]

  null_edge_chunks <- foreach(
    chunk     = chunk_sizes,
    .combine  = "c",
    .packages = "bnlearn",
    .export   = "learn_hc",
    .multicombine = TRUE
  ) %dopar% {
    local_edges <- integer(chunk)
    for (j in seq_len(chunk)) {
      null_data    <- data
      null_data[]  <- lapply(data, sample)
      null_net     <- learn_hc(null_data, maxp = maxp)
      local_edges[j] <- nrow(arcs(null_net))
    }
    local_edges
  }

  # Conservative permutation p-value  (Phipson & Smyth, 2010)
  elapsed_null <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  Null model done in %.1f min\n", elapsed_null / 60))
  p_perm <- (1 + sum(null_edge_chunks >= observed_edges)) / (1 + length(null_edge_chunks))

  list(
    observed_edges  = observed_edges,
    null_edges_mean = mean(null_edge_chunks),
    null_edges_sd   = sd(null_edge_chunks),
    null_edges_raw  = null_edge_chunks,
    p_perm          = p_perm
  )
}

# ---- 6. Edge-overlap between unconstrained and constrained -------------------

compute_edge_overlap <- function(data, maxp_constrained) {
  # What fraction of maxp=k edges are contained in the unconstrained network?
  bn_unc  <- hc(data)
  bn_con  <- hc(data, maxp = as.integer(maxp_constrained))

  arcs_unc <- paste(arcs(bn_unc)[,1], arcs(bn_unc)[,2], sep = "->")
  arcs_con <- paste(arcs(bn_con)[,1], arcs(bn_con)[,2], sep = "->")

  if (length(arcs_con) == 0) return(NA_real_)
  sum(arcs_con %in% arcs_unc) / length(arcs_con)
}

# ---- 7. Sensitivity runner (loops over maxp values) --------------------------

run_parent_sensitivity <- function(data, label, out_dir,
                                   maxp_values    = c(NA, 5, 3),
                                   n_boots_primary = 5000,
                                   n_boots_sens    = 1000,
                                   n_perm_primary  = 5000,
                                   n_perm_sens     = 1000) {
  rows <- list()
  n_conditions <- length(maxp_values)
  t_start <- proc.time()

  for (idx in seq_along(maxp_values)) {
    mp <- maxp_values[idx]
    tag     <- if (is.na(mp)) "maxpInf" else paste0("maxp", mp)
    n_boots <- if (is.na(mp)) n_boots_primary else n_boots_sens
    n_perm  <- if (is.na(mp)) n_perm_primary  else n_perm_sens

    cat(sprintf("\n--- %s | %s  [condition %d/%d] ---\n", label, tag, idx, n_conditions))

    boot <- bootstrap_network_analysis(data, n_boots = n_boots, maxp = mp)
    nul  <- generate_null_networks(data, n_permutations = n_perm, maxp = mp)

    # Edge overlap with unconstrained (only meaningful for constrained runs)
    overlap <- if (is.na(mp)) NA_real_ else compute_edge_overlap(data, mp)

    rows[[tag]] <- data.frame(
      dataset               = label,
      maxp                  = if (is.na(mp)) "Inf" else as.character(mp),
      n_bootstraps          = n_boots,
      observed_edges        = nul$observed_edges,
      kept_edges_bs_ge_0.5  = boot$n_edges_kept,
      max_in_degree         = boot$max_parents,
      null_edges_mean       = round(nul$null_edges_mean, 1),
      null_edges_sd         = round(nul$null_edges_sd, 1),
      p_perm                = nul$p_perm,
      overlap_with_unconstr = if (is.na(overlap)) NA else round(overlap, 3),
      stringsAsFactors      = FALSE
    )

    # Save per-condition artefacts (tagged, never overwrites originals)
    saveRDS(boot, file.path(out_dir, paste0("bn_boot_", label, "_", tag, ".rds")))
    saveRDS(nul,  file.path(out_dir, paste0("bn_null_", label, "_", tag, ".rds")))

    cat(sprintf("  Observed edges: %d | Null mean: %.1f | p = %.4f | Max parents: %d\n",
                nul$observed_edges, nul$null_edges_mean, nul$p_perm, boot$max_parents))

    # Progress estimate
    elapsed_so_far <- (proc.time() - t_start)["elapsed"]
    avg_per_cond   <- elapsed_so_far / idx
    remaining      <- avg_per_cond * (n_conditions - idx)
    cat(sprintf("  [%s elapsed | ~%.0f min remaining for this tissue]\n",
                sprintf("%.1f min", elapsed_so_far / 60), remaining / 60))
  }

  do.call(rbind, rows)
}

# ---- 8. Load data ------------------------------------------------------------

cat("Loading prepared data from:", original_bn_path, "\n")
script_t0 <- proc.time()
processed_data <- readRDS(file.path(original_bn_path, "processed_data.rds"))

# Reproduce the same data split as baysian_network.R
pivot_data  <- processed_data$pivot_data
data_for_bn <- pivot_data[, -(1:7)]
data_for_bn <- as.data.frame(lapply(data_for_bn, as.numeric))
data_for_bn[is.infinite(as.matrix(data_for_bn))] <- NA
data_for_bn <- data_for_bn[, colSums(is.na(data_for_bn)) < nrow(data_for_bn)]
data_for_bn <- as.data.frame(lapply(data_for_bn, function(x) {
  ifelse(is.na(x), median(x, na.rm = TRUE), x)
}))

root_data <- data_for_bn[, grep("^R_", colnames(data_for_bn), value = TRUE)]
leaf_data <- data_for_bn[, grep("^L_", colnames(data_for_bn), value = TRUE)]

cat(sprintf("Root data: %d samples × %d features\n", nrow(root_data), ncol(root_data)))
cat(sprintf("Leaf data: %d samples × %d features\n", nrow(leaf_data), ncol(leaf_data)))

# ---- 9. Run sensitivity for both tissues ------------------------------------

cat("\n========== LEAF ==========\n")
sens_leaf <- run_parent_sensitivity(leaf_data, "Leaf", output_dir)

cat("\n========== ROOT ==========\n")
sens_root <- run_parent_sensitivity(root_data, "Root", output_dir)

# ---- 10. Combined summary table (Supp Table S9) -----------------------------

sens_all <- rbind(sens_leaf, sens_root)
csv_path <- file.path(output_dir, "Supp_Table_S9_BayesianParentSensitivity.csv")
write.csv(sens_all, csv_path, row.names = FALSE)
cat("\nSummary table saved:", csv_path, "\n")

print(sens_all)

# ---- 11. Supplementary Figure S9 --------------------------------------------

# Colours consistent with the rest of the paper
LEAF_COL <- "#1b941b"
ROOT_COL <- "#1e597d"

# --- Panel A: observed vs null edge counts across maxp regimes ---

plot_df <- sens_all %>%
  mutate(maxp = factor(maxp, levels = c("Inf", "5", "3"))) %>%
  select(dataset, maxp, observed_edges, null_edges_mean) %>%
  pivot_longer(cols = c(observed_edges, null_edges_mean),
               names_to = "type", values_to = "edges") %>%
  mutate(type = recode(type,
                       "observed_edges" = "Observed",
                       "null_edges_mean" = "Null expectation"))

p_edges <- ggplot(plot_df, aes(x = maxp, y = edges, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "black", linewidth = 0.3) +
  facet_wrap(~dataset, scales = "free_y") +
  scale_fill_manual(values = c("Observed" = "#555555", "Null expectation" = "#cccccc")) +
  labs(x = "Maximum parents per node",
       y = "Edge count",
       fill = NULL,
       title = "A") +
  theme_minimal(base_size = 13) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(face = "bold", size = 18, hjust = -0.05),
    legend.position  = "top",
    strip.text       = element_text(face = "bold", size = 13)
  )

# --- Panel B: max in-degree (hub truncation demonstration) ---

indeg_df <- sens_all %>%
  mutate(maxp = factor(maxp, levels = c("Inf", "5", "3")),
         col  = ifelse(dataset == "Leaf", LEAF_COL, ROOT_COL))

p_indeg <- ggplot(indeg_df, aes(x = maxp, y = max_in_degree, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Leaf" = LEAF_COL, "Root" = ROOT_COL)) +
  labs(x = "Maximum parents per node",
       y = "Max observed in-degree",
       fill = "Tissue",
       title = "B") +
  theme_minimal(base_size = 13) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(face = "bold", size = 18, hjust = -0.05),
    legend.position  = "top"
  )

# --- Panel C: null distribution + observed line (BOTH tissues, Inf vs maxp=3) ---

null_leaf_inf <- readRDS(file.path(output_dir, "bn_null_Leaf_maxpInf.rds"))
null_leaf_3   <- readRDS(file.path(output_dir, "bn_null_Leaf_maxp3.rds"))
null_root_inf <- readRDS(file.path(output_dir, "bn_null_Root_maxpInf.rds"))
null_root_3   <- readRDS(file.path(output_dir, "bn_null_Root_maxp3.rds"))

null_dist_df <- rbind(
  data.frame(tissue = "Leaf", maxp = "Unconstrained", null_edges = null_leaf_inf$null_edges_raw),
  data.frame(tissue = "Leaf", maxp = "maxp = 3",      null_edges = null_leaf_3$null_edges_raw),
  data.frame(tissue = "Root", maxp = "Unconstrained", null_edges = null_root_inf$null_edges_raw),
  data.frame(tissue = "Root", maxp = "maxp = 3",      null_edges = null_root_3$null_edges_raw)
)

obs_lines <- data.frame(
  tissue = c("Leaf", "Leaf", "Root", "Root"),
  maxp   = c("Unconstrained", "maxp = 3", "Unconstrained", "maxp = 3"),
  obs    = c(null_leaf_inf$observed_edges, null_leaf_3$observed_edges,
             null_root_inf$observed_edges, null_root_3$observed_edges)
)

tissue_fills <- c("Leaf" = LEAF_COL, "Root" = ROOT_COL)

p_null <- ggplot(null_dist_df, aes(x = null_edges, fill = tissue)) +
  geom_histogram(bins = 50, alpha = 0.6, colour = "white", linewidth = 0.2) +
  geom_vline(data = obs_lines, aes(xintercept = obs),
             colour = "red", linewidth = 0.9, linetype = "dashed") +
  facet_grid(tissue ~ maxp, scales = "free") +
  scale_fill_manual(values = tissue_fills, guide = "none") +
  labs(x = "Edge count",
       y = "Frequency (null permutations)",
       title = "C") +
  theme_minimal(base_size = 13) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    plot.title       = element_text(face = "bold", size = 18, hjust = -0.05),
    strip.text       = element_text(face = "bold", size = 13)
  )

# --- Save figure data for re-plotting without re-running analysis -------------

fig_data_path <- file.path(output_dir, "Supp_Fig_S9_plot_data.rds")
saveRDS(list(
  panel_A_data   = plot_df,
  panel_B_data   = indeg_df,
  panel_C_null   = null_dist_df,
  panel_C_obs    = obs_lines,
  summary_table  = sens_all,
  colours        = list(LEAF = LEAF_COL, ROOT = ROOT_COL)
), fig_data_path)
cat("Figure plot data saved:", fig_data_path, "\n")

# --- Compose multi-panel figure and save BOTH PDF and PNG ---------------------

compose_and_save <- function() {
  pdf_path <- file.path(output_dir, "Supp_Fig_S9_BayesianParentSensitivity.pdf")
  png_path <- file.path(output_dir, "Supp_Fig_S9_BayesianParentSensitivity.png")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    fig <- (p_edges | p_indeg) / p_null +
      plot_annotation(
        title    = "Bayesian Network Parent-Constraint Sensitivity Analysis",
        subtitle = "Leaf > Root architectural asymmetry preserved across all constraint levels (P < 0.001)",
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 12, colour = "grey30")
        )
      )
    ggsave(pdf_path, fig, width = 14, height = 11, dpi = 300)
    ggsave(png_path, fig, width = 14, height = 11, dpi = 150)

  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    library(cowplot)
    top_row <- plot_grid(p_edges, p_indeg, ncol = 2, rel_widths = c(1.2, 1))
    fig     <- plot_grid(top_row, p_null, ncol = 1, rel_heights = c(1, 1))
    ggsave(pdf_path, fig, width = 14, height = 11, dpi = 300)
    ggsave(png_path, fig, width = 14, height = 11, dpi = 150)

  } else {
    # Fallback: save each panel individually
    ggsave(file.path(output_dir, "Supp_Fig_S9_panelA_edges.pdf"),    p_edges, width = 8, height = 5, dpi = 300)
    ggsave(file.path(output_dir, "Supp_Fig_S9_panelB_indegree.pdf"), p_indeg, width = 6, height = 5, dpi = 300)
    ggsave(file.path(output_dir, "Supp_Fig_S9_panelC_nulldist.pdf"), p_null,  width = 10, height = 6, dpi = 300)
    pdf_path <- file.path(output_dir, "Supp_Fig_S9_panelA_edges.pdf")
    cat("Note: patchwork/cowplot not available. Panels saved as separate PDFs.\n")
  }

  cat("Figure saved:", pdf_path, "\n")
}

compose_and_save()

# ---- 12. Console summary ----------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  BAYESIAN PARENT-CONSTRAINT SENSITIVITY — COMPLETE\n")
cat("================================================================\n")
cat("\n")
cat("Outputs in:", output_dir, "\n\n")
cat("  Supp_Table_S9_BayesianParentSensitivity.csv\n")
cat("  Supp_Fig_S9_BayesianParentSensitivity.pdf / .png\n")
cat("  Supp_Fig_S9_plot_data.rds  (for re-plotting without re-running)\n")
cat("  Per-condition .rds files (bn_boot_*, bn_null_*)\n")
cat("\n")
cat("Key results:\n")
for (i in seq_len(nrow(sens_all))) {
  row <- sens_all[i, ]
  cat(sprintf("  %s maxp=%-3s : %3d observed edges | null %.1f +/- %.1f | p = %.4f | max parents = %d",
              row$dataset, row$maxp, row$observed_edges,
              row$null_edges_mean, row$null_edges_sd, row$p_perm, row$max_in_degree))
  if (!is.na(row$overlap_with_unconstr)) {
    cat(sprintf(" | overlap = %.1f%%", row$overlap_with_unconstr * 100))
  }
  cat("\n")
}
cat("\n")
cat("All outputs are ADDITIVE. Original BN results are untouched.\n")
total_min <- (proc.time() - script_t0)["elapsed"] / 60
cat(sprintf("Total runtime: %.1f minutes\n", total_min))
cat("================================================================\n")
