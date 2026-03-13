#!/usr/bin/env python3
"""
bootstrap_stability.py v2 — Cluster bootstrap + jackknife network metric
stability for TPC-2025-1173 revision.

Addresses:
  R1-Meth-1  "Are 7 biological replicates sufficient?"

Approach:
  (i)  Leave-one-replicate-out jackknife: drop all rows for replicate i,
       rebuild network → shows per-replicate influence on metrics
  (ii) Cluster bootstrap (200×): draw 7 replicate IDs with replacement,
       keep ALL rows from selected replicates, rebuild network
       → provides 95% CIs for density, transitivity, modularity
  Key test: Do leaf and root CIs remain non-overlapping for density and
  modularity? If yes, the architectural asymmetry is robust to sampling
  variation at n=7 biological replicates.

Design decisions:
  - CLUSTER resampling by biological replicate (column "Replication"),
    not by row. Each genotype has 84 rows = 7 replicates × 2 treatments
    × 6 timepoints. The reviewer question is about replicate sufficiency,
    so the resampling unit must be the replicate.
  - Full BH-FDR in every iteration (HANDOVER principle #2)
  - Metrics on both full graph and LCC (HANDOVER principle #4)
  - weight=|ρ| for modularity, rho_signed preserved (v2 fix #4)
  - Independent RNG seeds per tissue×genotype group
  - python-louvain is a hard dependency (we cite modularity CIs)

Outputs:
  Supp Table S6  — Per-iteration metrics (baseline + jackknife + bootstrap)
  Supp Table S6b — Summary: baseline, median, 95% CI, jackknife CV, overlap test
  Supp Fig  S6   — Bootstrap distributions with 95% CIs for 3 key metrics

Usage:
  python bootstrap_stability.py --data-root "<data_root>" --out "output/revision_pack"

  Optional:
    --vip-combined-rel  (default: output\\results\\spearman\\VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv)
    --rho       (default: 0.7)
    --fdr       (default: 0.05)
    --n-boot    (default: 200)
    --seed      (default: 42)
    --n-workers (default: 0 = auto)
    --rep-col   (default: Replication)
    --no-plot
"""

import argparse
import os
import sys
import time
import warnings
from pathlib import Path
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import rankdata, t as t_dist

try:
    import community as community_louvain
except ImportError:
    sys.exit(
        "[ERROR] python-louvain is required (modularity CIs are a key output).\n"
        "Install with: pip install python-louvain"
    )

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib not installed. Figures will be skipped.")

warnings.filterwarnings("ignore", category=RuntimeWarning)

# Limit BLAS threading when using process-level parallelism
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")


# ============================================================================
# 1. CORE FUNCTIONS (shared with network_robustness.py v2)
# ============================================================================

def bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out


def build_network(data: np.ndarray, col_names: list,
                  rho_thr: float, fdr_alpha: float) -> nx.Graph:
    """Build Spearman correlation network with full BH-FDR.

    Optimised: Fast rank-Pearson via np.corrcoef, 1D upper-triangle p-values,
    bulk edge insertion. Mathematically identical to scipy.stats.spearmanr
    but ~100× faster for large feature matrices.

    Edge attributes: weight=|ρ| (for Louvain), rho_signed=ρ (biology).
    """
    n_samples, n_features = data.shape
    G = nx.Graph()
    G.add_nodes_from(col_names)

    if n_features < 2 or n_samples < 3:
        return G

    # 1. Rank data column-wise → Spearman = Pearson(ranks)
    ranked_data = rankdata(data, method='average', axis=0)

    # 2. Pearson correlation on ranks (uses LAPACK, much faster than scipy spearmanr)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        corr = np.corrcoef(ranked_data, rowvar=False)

    # 3. Extract upper triangle immediately (work in 1D from here)
    iu = np.triu_indices(n_features, k=1)
    rho_vals = corr[iu]

    # 4. Handle NaNs from zero-variance columns → force p=1.0
    valid_mask = ~np.isnan(rho_vals)

    # 5. Vectorised p-value calculation (t-distribution, only on valid pairs)
    pvals = np.ones_like(rho_vals)  # NaN/invalid → p=1.0 (excluded by BH)
    r_sq = np.clip(rho_vals[valid_mask] ** 2, 0, 1.0)
    t_stat = rho_vals[valid_mask] * np.sqrt(
        (n_samples - 2) / ((1.0 - r_sq) + 1e-15))
    pvals[valid_mask] = 2 * t_dist.sf(np.abs(t_stat), df=n_samples - 2)

    # 6. BH-FDR and thresholding
    qvals = bh_qvalues(pvals)
    keep_mask = (qvals <= fdr_alpha) & (np.abs(rho_vals) > rho_thr) & valid_mask
    keep_idx = np.where(keep_mask)[0]

    # 7. Bulk edge insertion (avoids per-edge Python interpreter overhead)
    edges_to_add = [
        (col_names[iu[0][idx]], col_names[iu[1][idx]],
         {"weight": float(abs(rho_vals[idx])), "rho_signed": float(rho_vals[idx])})
        for idx in keep_idx
    ]
    G.add_edges_from(edges_to_add)

    return G


def lcc_subgraph(G: nx.Graph):
    """Extract largest connected component."""
    if G.number_of_edges() == 0:
        return None
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    if not comps:
        return None
    return G.subgraph(comps[0]).copy()


def modularity_louvain(G: nx.Graph) -> float:
    """Louvain modularity using weight=|ρ| (non-negative)."""
    if G.number_of_edges() == 0 or G.number_of_nodes() < 2:
        return float("nan")
    try:
        part = community_louvain.best_partition(G, weight="weight", random_state=42)
        return float(community_louvain.modularity(part, G, weight="weight"))
    except Exception:
        return float("nan")


def safe_mean_path_length(G: nx.Graph) -> float:
    if G.number_of_edges() == 0 or G.number_of_nodes() < 2:
        return float("nan")
    try:
        return float(nx.average_shortest_path_length(G))
    except Exception:
        return float("nan")


def metrics_bundle(G: nx.Graph) -> dict:
    """Compute metrics on both full graph and LCC."""
    n_all = int(G.number_of_nodes())
    n_iso = int(sum(1 for _, d in G.degree() if d == 0))
    e_all = int(G.number_of_edges())
    dens_all = float(nx.density(G)) if n_all > 1 else float("nan")

    H = lcc_subgraph(G)
    if H is None:
        return {
            "n_nodes_all": n_all, "n_isolates": n_iso, "n_edges": e_all,
            "density_all": dens_all,
            "modularity_all": float("nan"),
            "lcc_nodes": 0, "density_lcc": float("nan"),
            "transitivity_lcc": float("nan"), "modularity_lcc": float("nan"),
            "path_length_lcc": float("nan"),
        }

    return {
        "n_nodes_all": n_all, "n_isolates": n_iso, "n_edges": e_all,
        "density_all": dens_all,
        "modularity_all": modularity_louvain(G),  # Full graph (matching manuscript)
        "lcc_nodes": int(H.number_of_nodes()),
        "density_lcc": float(nx.density(H)) if H.number_of_nodes() > 1 else float("nan"),
        "transitivity_lcc": float(nx.transitivity(H)) if H.number_of_nodes() > 2 else float("nan"),
        "modularity_lcc": modularity_louvain(H),
        "path_length_lcc": float("nan"),  # BYPASSED: O(V×E) bottleneck, not in reviewer test
    }


# ============================================================================
# 2. POOL INITIALIZER + CLUSTER RESAMPLING WORKER
# ============================================================================

# Global variable holding data slices per worker process.
# Populated once via Pool initializer — avoids pickling data per task.
_WORKER_SLICES = {}


def _init_worker(shared_slices: dict):
    """Pool initializer: runs once per worker process."""
    global _WORKER_SLICES
    _WORKER_SLICES = shared_slices


def _worker(task: dict) -> dict:
    """Process one baseline / jackknife / bootstrap iteration.

    Cluster resampling: the unit is the biological replicate, not the row.
    Big data arrays are read from _WORKER_SLICES (set once per process),
    not from the task dict — eliminates IPC pickle overhead.
    """
    tissue = task["tissue"]
    genotype = task["genotype"]

    # Pull data from worker-global memory (set by _init_worker)
    s = _WORKER_SLICES[(tissue, genotype)]
    data = s["data"]
    col_names = s["col_names"]
    rep_ids = s["rep_ids"]
    unique_reps = s["unique_reps"]

    rho_thr = task["rho_thr"]
    fdr_alpha = task["fdr_alpha"]
    resample = task["resample"]
    iter_idx = task["iter_idx"]

    n_reps = len(unique_reps)

    if resample == "baseline":
        X = data

    elif resample == "jackknife":
        # Leave out all rows from replicate at index iter_idx
        drop_rep = unique_reps[iter_idx]
        keep_mask = rep_ids != drop_rep
        X = data[keep_mask]

    elif resample == "bootstrap":
        # Draw n_reps replicate IDs with replacement, keep all their rows
        rng = np.random.default_rng(task["seed"])
        sampled_reps = rng.choice(unique_reps, size=n_reps, replace=True)
        row_indices = []
        for rep in sampled_reps:
            row_indices.append(np.where(rep_ids == rep)[0])
        row_indices = np.concatenate(row_indices)
        X = data[row_indices]
    else:
        raise ValueError(f"Unknown resample type: {resample}")

    G = build_network(X, col_names, rho_thr, fdr_alpha)
    m = metrics_bundle(G)
    m.update({
        "tissue": tissue,
        "genotype": genotype,
        "resample": resample,
        "iter": iter_idx,
        "n_rows_used": X.shape[0],
        "n_reps_used": n_reps if resample != "jackknife" else n_reps - 1,
    })
    return m


# ============================================================================
# 3. FIGURE GENERATION
# ============================================================================

def make_figure(df: pd.DataFrame, out_dir: Path):
    """Supp Fig S6: Bootstrap metric distributions with 95% CIs."""
    if not HAS_MPL:
        print("[SKIP] matplotlib not available.")
        return

    boot = df[df["resample"] == "bootstrap"].copy()
    if boot.empty:
        print("[SKIP] No bootstrap data for figure.")
        return

    metrics = [
        ("density_all", "Density (full graph)"),
        ("modularity_all", "Modularity (full graph)"),
        ("transitivity_lcc", "Transitivity (LCC)"),
    ]

    color_map = {
        ("Leaf", "G1"):  "#2ca02c",
        ("Leaf", "G2"):  "#98df8a",
        ("Root", "G1"):  "#d62728",
        ("Root", "G2"):  "#ff9896",
    }
    groups = [("Leaf", "G1"), ("Leaf", "G2"), ("Root", "G1"), ("Root", "G2")]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        "Supplementary Figure S6: Cluster-bootstrap network metric stability\n"
        "(n = 7 biological replicates, 200 resamples, full BH-FDR)",
        fontsize=11, fontweight="bold", y=1.04)

    for ax, (metric_col, metric_label) in zip(axes, metrics):
        positions = []
        box_data = []
        colors = []
        labels = []
        baselines = []

        for i, (tissue, genotype) in enumerate(groups):
            vals = boot[(boot["tissue"] == tissue) & (boot["genotype"] == genotype)][metric_col].dropna()
            if len(vals) == 0:
                continue
            positions.append(i)
            box_data.append(vals.values)
            colors.append(color_map[(tissue, genotype)])
            labels.append(f"{tissue}\n{genotype}")

            # Baseline value
            base = df[(df["tissue"] == tissue) & (df["genotype"] == genotype) &
                      (df["resample"] == "baseline")][metric_col]
            if not base.empty:
                baselines.append((i, float(base.iloc[0])))

        if box_data:
            bp = ax.boxplot(box_data, positions=positions, widths=0.5,
                           patch_artist=True, showfliers=False, zorder=3)
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.6)
            for element in ["whiskers", "caps"]:
                for line in bp[element]:
                    line.set_color("gray")
                    line.set_linewidth(1)
            for line in bp["medians"]:
                line.set_color("black")
                line.set_linewidth(1.5)

        # Mark baseline
        for pos, val in baselines:
            ax.plot(pos, val, marker="*", color="black", markersize=10, zorder=5)

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel(metric_label, fontsize=10)
        ax.set_title(metric_label, fontsize=10)
        ax.grid(axis="y", alpha=0.3)

    # Legend for baseline marker
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker="*", color="black", linestyle="None",
                              markersize=10, label="Baseline (all 7 replicates)")]
    axes[2].legend(handles=legend_elements, loc="upper right", fontsize=8)

    fig.tight_layout()
    out_path = out_dir / "Supp_Fig_S6_WheatNetworkMetricStability.pdf"
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  [OK] {out_path}")


# ============================================================================
# 4. SUMMARY TABLE
# ============================================================================

def compute_summary(df: pd.DataFrame) -> tuple:
    """Compute summary statistics from bootstrap results.

    Returns:
        summary_df: per-metric CI and stability stats
        paired_df:  paired difference CI (leaf − root, per bootstrap iteration)
    """
    metrics = ["density_all", "modularity_all", "density_lcc", "transitivity_lcc",
               "modularity_lcc", "path_length_lcc", "n_edges", "n_isolates"]

    rows = []
    for tissue in ["Leaf", "Root"]:
        for genotype in ["G1", "G2"]:
            mask_base = (df["tissue"] == tissue) & (df["genotype"] == genotype)
            base_row = df[mask_base & (df["resample"] == "baseline")]
            boot_rows = df[mask_base & (df["resample"] == "bootstrap")]
            jack_rows = df[mask_base & (df["resample"] == "jackknife")]

            for m in metrics:
                if m not in df.columns:
                    continue
                base_val = float(base_row[m].iloc[0]) if not base_row.empty else float("nan")
                boot_vals = boot_rows[m].dropna().values

                if len(boot_vals) > 0:
                    ci_lo = float(np.percentile(boot_vals, 2.5))
                    ci_hi = float(np.percentile(boot_vals, 97.5))
                    boot_med = float(np.median(boot_vals))
                    boot_std = float(np.std(boot_vals))
                    ci_width = ci_hi - ci_lo
                else:
                    ci_lo = ci_hi = boot_med = boot_std = ci_width = float("nan")

                jack_vals = jack_rows[m].dropna().values
                if len(jack_vals) > 1 and np.mean(jack_vals) != 0:
                    jack_cv = float(np.std(jack_vals) / abs(np.mean(jack_vals)))
                else:
                    jack_cv = float("nan")

                rows.append({
                    "tissue": tissue, "genotype": genotype, "metric": m,
                    "baseline": base_val,
                    "bootstrap_median": boot_med,
                    "bootstrap_ci_lo": ci_lo,
                    "bootstrap_ci_hi": ci_hi,
                    "bootstrap_ci_width": ci_width,
                    "bootstrap_std": boot_std,
                    "jackknife_cv": jack_cv,
                })

    summary = pd.DataFrame(rows)

    # ── PAIRED DIFFERENCE TEST ──
    # Because leaf and root bootstrap iterations share the same seed per
    # genotype, iteration b samples the same replicate IDs for both tissues.
    # This makes the difference Δ(b) = leaf(b) − root(b) a proper paired
    # statistic. Testing whether the CI of Δ excludes zero is more powerful
    # than checking individual CI overlap.
    paired_rows = []
    test_metrics = [
        ("density_all", "leaf_denser", "positive"),     # expect Leaf > Root → Δ > 0
        ("modularity_all", "root_more_modular", "negative"),  # expect Root > Leaf → Δ < 0
    ]

    for genotype in ["G1", "G2"]:
        for m, test_name, expected_sign in test_metrics:
            if m not in df.columns:
                continue
            leaf_boot = (df[(df["genotype"] == genotype) & (df["tissue"] == "Leaf") &
                           (df["resample"] == "bootstrap")]
                        .sort_values("iter").reset_index(drop=True))
            root_boot = (df[(df["genotype"] == genotype) & (df["tissue"] == "Root") &
                           (df["resample"] == "bootstrap")]
                        .sort_values("iter").reset_index(drop=True))

            if len(leaf_boot) == 0 or len(root_boot) == 0:
                continue

            # Paired difference
            diff = leaf_boot[m].values - root_boot[m].values
            diff = diff[~np.isnan(diff)]

            if len(diff) == 0:
                continue

            ci_lo = float(np.percentile(diff, 2.5))
            ci_hi = float(np.percentile(diff, 97.5))
            mean_diff = float(np.mean(diff))

            # Does CI exclude zero?
            if expected_sign == "positive":
                significant = ci_lo > 0  # Entire CI above zero
            else:
                significant = ci_hi < 0  # Entire CI below zero

            # Baseline difference
            base_leaf = df[(df["genotype"] == genotype) & (df["tissue"] == "Leaf") &
                          (df["resample"] == "baseline")]
            base_root = df[(df["genotype"] == genotype) & (df["tissue"] == "Root") &
                          (df["resample"] == "baseline")]
            if not base_leaf.empty and not base_root.empty:
                base_diff = float(base_leaf[m].iloc[0]) - float(base_root[m].iloc[0])
            else:
                base_diff = float("nan")

            paired_rows.append({
                "genotype": genotype,
                "metric": m,
                "test": test_name,
                "expected_sign": expected_sign,
                "baseline_diff": base_diff,
                "mean_boot_diff": mean_diff,
                "ci_lo": ci_lo,
                "ci_hi": ci_hi,
                "ci_excludes_zero": significant,
                "n_boot_pairs": len(diff),
            })

    return summary, pd.DataFrame(paired_rows)


# ============================================================================
# 5. MAIN
# ============================================================================

def main():
    t0 = time.time()

    ap = argparse.ArgumentParser(
        description="Cluster bootstrap + jackknife network metric stability (TPC-2025-1173)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--leaf-data-rel", default=os.path.join("data", "data", "n_p_l.csv"))
    ap.add_argument("--root-data-rel", default=os.path.join("data", "data", "n_p_r.csv"))
    ap.add_argument("--vip-combined-rel",
                    default=os.path.join("output", "results", "spearman",
                                         "VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv"))
    ap.add_argument("--rho", type=float, default=0.7)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-workers", type=int, default=0,
                    help="Parallel workers (0=auto, 1=sequential)")
    ap.add_argument("--rep-col", default="Replication",
                    help="Column name for biological replicate ID")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    data_root = Path(args.data_root)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.n_workers <= 0:
        n_workers = max(1, cpu_count() - 2)
    else:
        n_workers = args.n_workers

    print("=" * 70)
    print("CLUSTER BOOTSTRAP / JACKKNIFE STABILITY — TPC-2025-1173 Revision")
    print("=" * 70)
    print(f"  Data root:       {data_root}")
    print(f"  Output dir:      {out_dir}")
    print(f"  |rho|:           {args.rho}")
    print(f"  FDR alpha:       {args.fdr}")
    print(f"  Bootstrap:       {args.n_boot} iterations (cluster by {args.rep_col})")
    print(f"  Seed:            {args.seed}")
    print(f"  Workers:         {n_workers}")
    print(f"  BH-FDR:          Full (in every iteration)")
    print(f"  Resampling unit: Biological replicate (not row)")
    print()

    # ── Load data ──
    print("[1/4] Loading data...")
    leaf_path = data_root / args.leaf_data_rel
    root_path = data_root / args.root_data_rel
    vip_path = data_root / args.vip_combined_rel

    for p, name in [(leaf_path, "Leaf"), (root_path, "Root"), (vip_path, "VIP")]:
        if not p.exists():
            sys.exit(f"[ERROR] {name} file not found: {p}")

    leaf = pd.read_csv(leaf_path)
    root = pd.read_csv(root_path)
    vip = pd.read_csv(vip_path)
    vip = vip.rename(columns={vip.columns[0]: "Metabolite"})
    vip_set = set(vip["Metabolite"].astype(str))

    # Verify replicate column exists
    for tissue_name, df in [("Leaf", leaf), ("Root", root)]:
        if args.rep_col not in df.columns:
            sys.exit(f"[ERROR] Replicate column '{args.rep_col}' not found in {tissue_name} data.\n"
                     f"  Available columns: {list(df.columns[:15])}")

    print(f"  Leaf: {leaf.shape[0]} rows x {leaf.shape[1] - 10} features")
    print(f"  Root: {root.shape[0]} rows x {root.shape[1] - 10} features")
    print(f"  VIP features (>1.0): {len(vip_set)}")

    # ── Prepare data slices per (tissue, genotype) ──
    # Filter to Treatment=1 (stressed) only — matching manuscript network construction.
    # Original spearman_network.py filters treatment before building networks.
    slices = {}
    for tissue_name, df in [("Leaf", leaf), ("Root", root)]:
        X_all = df.iloc[:, 10:]
        meta = df[["Genotype", "Treatment", args.rep_col]].copy()
        cols = [c for c in X_all.columns if c in vip_set]

        for genotype in ["G1", "G2"]:
            mask = (meta["Genotype"].astype(str) == genotype) & (meta["Treatment"] == 1)
            X = X_all.loc[mask, cols].copy()
            rep_ids = meta.loc[mask, args.rep_col].values

            X = X.dropna(axis=1, how="all")
            # Keep rep_ids aligned after any potential row drops (shouldn't happen here)

            unique_reps = sorted(set(rep_ids))
            n_reps = len(unique_reps)

            if X.shape[1] < 3 or n_reps < 3:
                print(f"  [WARN] {tissue_name} {genotype}: {n_reps} replicates, "
                      f"{X.shape[1]} features — skipping")
                continue

            slices[(tissue_name, genotype)] = {
                "data": X.values,
                "col_names": list(X.columns),
                "rep_ids": rep_ids,
                "unique_reps": unique_reps,
                "n_rows": X.shape[0],
                "n_features": X.shape[1],
                "n_reps": n_reps,
            }
            rows_per_rep = X.shape[0] // n_reps
            print(f"  {tissue_name} {genotype}: {n_reps} replicates x "
                  f"{rows_per_rep} rows/rep = {X.shape[0]} rows, "
                  f"{X.shape[1]} features")
    print()

    # ── Build task list (PAIRED seeds per genotype) ──
    print("[2/4] Building task list...")

    # Master RNG for generating per-genotype seed streams.
    # Within a genotype, leaf and root bootstrap iter b use the SAME seed
    # → same replicates sampled → enables paired difference test.
    master_rng = np.random.default_rng(args.seed)
    genotype_seeds = {}
    for genotype in ["G1", "G2"]:
        genotype_seeds[genotype] = master_rng.integers(0, 2**31, size=args.n_boot)

    tasks = []

    for (tissue, genotype), s in slices.items():
        common = {
            "rho_thr": args.rho,
            "fdr_alpha": args.fdr,
            "tissue": tissue,
            "genotype": genotype,
        }

        # Baseline
        tasks.append({**common, "resample": "baseline", "iter_idx": 0, "seed": 0})

        # Jackknife (leave-one-replicate-out)
        for i in range(s["n_reps"]):
            tasks.append({**common, "resample": "jackknife", "iter_idx": i, "seed": 0})

        # Bootstrap — uses genotype-level seeds (paired across tissues)
        for b in range(args.n_boot):
            tasks.append({**common, "resample": "bootstrap",
                         "iter_idx": b + 1, "seed": int(genotype_seeds[genotype][b])})

    n_baseline = sum(1 for t in tasks if t["resample"] == "baseline")
    n_jackknife = sum(1 for t in tasks if t["resample"] == "jackknife")
    n_bootstrap = sum(1 for t in tasks if t["resample"] == "bootstrap")
    print(f"  Baseline:       {n_baseline}")
    print(f"  Jackknife:      {n_jackknife} (leave-one-replicate-out, {n_jackknife // max(n_baseline,1)} per group)")
    print(f"  Bootstrap:      {n_bootstrap} (cluster by replicate)")
    print(f"  Total:          {len(tasks)} network builds")
    print()

    # ── Run ──
    print(f"[3/4] Running ({n_workers} workers)...")

    if n_workers == 1:
        # Sequential — manually initialise worker global state
        _init_worker(slices)
        results = []
        for i, task in enumerate(tasks):
            results.append(_worker(task))
            if (i + 1) % 20 == 0 or (i + 1) == len(tasks):
                print(f"  [{i+1}/{len(tasks)}] done")
    else:
        # Parallel — Pool initializer passes slices once per worker process
        with Pool(processes=n_workers,
                  initializer=_init_worker,
                  initargs=(slices,)) as pool:
            results = []
            for i, result in enumerate(pool.imap_unordered(_worker, tasks, chunksize=10)):
                results.append(result)
                if (i + 1) % 50 == 0 or (i + 1) == len(tasks):
                    print(f"  [{i+1}/{len(tasks)}] done")

    df = pd.DataFrame(results)

    # Save full iteration-level table
    out_s6 = out_dir / "Supp_Table_S6_WheatNetworkMetricStability.csv"
    df.to_csv(out_s6, index=False)
    print(f"  [OK] {out_s6}")
    print(f"       {len(df)} rows: {n_baseline} baseline + {n_jackknife} jackknife + {n_bootstrap} bootstrap")
    print()

    # ── Summary ──
    summary_df, overlap_df = compute_summary(df)

    out_summary = out_dir / "Supp_Table_S6b_WheatBootstrapSummary.csv"
    summary_df.to_csv(out_summary, index=False)
    print(f"  [OK] {out_summary}")

    print()
    print("  CLUSTER BOOTSTRAP 95% CIs (resampling unit = biological replicate):")
    for m in ["density_all", "modularity_all"]:
        print(f"    {m}:")
        for tissue in ["Leaf", "Root"]:
            for genotype in ["G1", "G2"]:
                row = summary_df[(summary_df["tissue"] == tissue) &
                                (summary_df["genotype"] == genotype) &
                                (summary_df["metric"] == m)]
                if not row.empty:
                    r = row.iloc[0]
                    print(f"      {tissue} {genotype}: {r['baseline']:.4f} "
                          f"[{r['bootstrap_ci_lo']:.4f}, {r['bootstrap_ci_hi']:.4f}] "
                          f"(jackknife CV={r['jackknife_cv']:.4f})")

    print()
    print("  PAIRED DIFFERENCE TEST (Leaf - Root, per bootstrap iteration):")
    print("  (Same replicates sampled for both tissues -> proper paired statistic)")
    if not overlap_df.empty:
        for _, row in overlap_df.iterrows():
            status = "SIGNIFICANT" if row["ci_excludes_zero"] else "NOT SIGNIFICANT"
            print(f"    {row['genotype']} {row['metric']} ({row['test']}): "
                  f"delta = {row['mean_boot_diff']:.4f} "
                  f"[{row['ci_lo']:.4f}, {row['ci_hi']:.4f}] -> {status}")
    else:
        print("    [WARN] Could not compute paired difference test.")
    print()

    # ── Figure ──
    if not args.no_plot:
        print("[4/4] Generating Supp Fig S6...")
        make_figure(df, out_dir)
    else:
        print("[4/4] Figure generation skipped (--no-plot).")

    elapsed = time.time() - t0
    print()
    print("=" * 70)
    print(f"DONE in {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"All outputs in: {out_dir}")
    print()
    print("Output files:")
    print(f"  Supp Table S6:  {out_s6.name}")
    print(f"  Supp Table S6b: {out_summary.name}")
    if not args.no_plot:
        print(f"  Supp Fig  S6:   Supp_Fig_S6_WheatNetworkMetricStability.pdf")
    print("=" * 70)


if __name__ == "__main__":
    main()