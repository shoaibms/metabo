#!/usr/bin/env python3
"""
network_robustness.py v2 — Consolidated network robustness analysis for TPC-2025-1173 revision.

Addresses:
  R2-1  Threshold sensitivity figure
  R2-4  VIP > 0.8 sensitivity
  R1-Meth-6  Singletons / disconnected components audit

Consolidates checklist sections: 2.1 (sweep) + 2.2 (plotting) + 2.7 (singleton audit)

Outputs:
  Supp Table S1  — Network metrics across VIP × |ρ| thresholds
  Supp Table S2  — Hub Jaccard overlap across thresholds
  Supp Table S7  — Singleton-missingness cross-reference audit
  Supp Fig  S1   — Threshold sensitivity (density, modularity, hub Jaccard vs |ρ|)
  Supp Fig  S2   — VIP sensitivity (metrics at VIP 0.8 vs 1.0)

Architecture:
  "Compute once, filter many" — Spearman correlation + p-value matrices are computed
  once per (VIP threshold × tissue × genotype). The ρ threshold sweep then just
  filters the precomputed matrices (a cheap mask operation), eliminating redundant
  correlation recomputation and avoiding multiprocessing pickle overhead on Windows.

Usage:
  python network_robustness.py --data-root "<data_root>" --out "output/revision_pack"

  Optional:
    --leaf-data-rel   (default: data\\data\\n_p_l.csv)
    --root-data-rel   (default: data\\data\\n_p_r.csv)
    --vip-leaf-rel    (default: output\\results\\updated_PLS\\vip_scores_all_leaf.csv)
    --vip-root-rel    (default: output\\results\\updated_PLS\\vip_scores_all_root.csv)
    --vip-combined-rel (default: output\\results\\spearman\\VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv)
    --vip-thresholds  (default: 1.0,0.8)
    --rho-thresholds  (default: 0.6,0.65,0.7,0.75,0.8)
    --fdr-alpha       (default: 0.05)
    --topk            (default: 20)
    --no-plot         (skip figure generation)
    --missing-s3b     (path to Supp_Table_S3b_MissingnessPerFeature.csv; if omitted,
                       missingness is computed on-the-fly from raw data)
"""

import argparse
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr

try:
    import community as community_louvain  # python-louvain
except ImportError:
    community_louvain = None
    print("[WARN] python-louvain not installed. Modularity will be NaN. "
          "Install with: pip install python-louvain")

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend for server/headless
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib not installed. Figures will be skipped.")

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ============================================================================
# 1. CORE NETWORK FUNCTIONS (adapted from spearman_network.py, FILE 30)
# ============================================================================

def bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values. Matches manuscript's FDR method (R2-6b)."""
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


def spearman_corr_p(df: pd.DataFrame):
    """Spearman correlation + p-value matrices across columns.

    Handles edge cases:
      - 0 or 1 feature: returns identity/zero matrices
      - 2 features: scipy returns scalar, we reshape to 2x2
      - 3+ features: standard matrix output
    """
    n_features = df.shape[1]
    cols = df.columns

    if n_features < 2:
        corr_df = pd.DataFrame(np.eye(n_features), index=cols, columns=cols)
        p_df = pd.DataFrame(np.zeros((n_features, n_features)), index=cols, columns=cols)
        return corr_df, p_df

    if n_features == 2:
        rho, pval = spearmanr(df.values[:, 0], df.values[:, 1])
        corr = np.array([[1.0, rho], [rho, 1.0]])
        ps = np.array([[0.0, pval], [pval, 0.0]])
        return pd.DataFrame(corr, index=cols, columns=cols), pd.DataFrame(ps, index=cols, columns=cols)

    corr, p = spearmanr(df.values, axis=0)
    return pd.DataFrame(corr, index=cols, columns=cols), pd.DataFrame(p, index=cols, columns=cols)


def build_graph_from_corr(corr_df: pd.DataFrame, p_df: pd.DataFrame,
                          rho_thr: float, fdr_alpha: float) -> nx.Graph:
    """Build network with BH-FDR correction + |rho| threshold.

    Full BH-FDR in every analysis per design principle #2.

    Edge attributes:
      weight     = |rho|  (non-negative; used by Louvain modularity)
      rho_signed = rho    (preserves biology: positive vs negative correlation)
    """
    cols = corr_df.columns
    n = len(cols)

    G = nx.Graph()
    G.add_nodes_from(cols)  # all VIP features are nodes (for singleton audit)

    if n < 2:
        return G

    iu = np.triu_indices(n, k=1)
    pvals = p_df.values[iu]
    qvals = bh_qvalues(pvals)

    sig = (qvals <= fdr_alpha)
    rho = corr_df.values[iu]
    keep = sig & (np.abs(rho) > rho_thr)

    edges_i = iu[0][keep]
    edges_j = iu[1][keep]
    weights = rho[keep]

    for i, j, w in zip(edges_i, edges_j, weights):
        G.add_edge(cols[i], cols[j], weight=float(abs(w)), rho_signed=float(w))
    return G


def modularity_louvain(G: nx.Graph) -> float:
    """Louvain modularity. Relies on weight=|rho| (non-negative) edge attribute."""
    if community_louvain is None or G.number_of_edges() == 0:
        return float("nan")
    try:
        part = community_louvain.best_partition(G, weight="weight", random_state=42)
        return float(community_louvain.modularity(part, G, weight="weight"))
    except Exception:
        return float("nan")


def safe_mean_path_length(G: nx.Graph) -> float:
    """Average shortest path length on the LCC."""
    if G.number_of_nodes() < 2 or G.number_of_edges() == 0:
        return float("nan")
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    if not comps:
        return float("nan")
    H = G.subgraph(comps[0]).copy()
    if H.number_of_nodes() < 2:
        return float("nan")
    try:
        return float(nx.average_shortest_path_length(H))
    except Exception:
        return float("nan")


def topk_hubs(G: nx.Graph, k: int):
    """Return top-k hub names by degree, excluding isolates (degree 0).

    Hubs are defined as the highest-degree *connected* nodes. Including
    degree-0 nodes would pollute hub sets and make Jaccard overlap
    meaningless.
    """
    deg = [(n, d) for n, d in G.degree() if d > 0]
    if not deg:
        return []
    deg.sort(key=lambda x: x[1], reverse=True)
    return [n for n, _ in deg[:k]]


def jaccard(a, b):
    """Jaccard similarity between two sets."""
    A, B = set(a), set(b)
    if len(A | B) == 0:
        return float("nan")
    return float(len(A & B) / len(A | B))


def component_stats(G: nx.Graph) -> dict:
    """Compute component-level statistics for singleton audit (R1-Meth-6)."""
    n_nodes = G.number_of_nodes()
    if n_nodes == 0:
        return {
            "n_components": 0, "n_singletons": 0, "n_components_2to5": 0,
            "lcc_nodes": 0, "lcc_frac": float("nan"),
        }

    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    comp_sizes = [len(c) for c in comps]
    lcc = comp_sizes[0] if comp_sizes else 0
    n_singletons = sum(1 for s in comp_sizes if s == 1)
    n_small = sum(1 for s in comp_sizes if 2 <= s <= 5)

    return {
        "n_components": len(comps),
        "n_singletons": n_singletons,
        "n_components_2to5": n_small,
        "lcc_nodes": int(lcc),
        "lcc_frac": float(lcc / n_nodes) if n_nodes else float("nan"),
    }


# ============================================================================
# 2. COMPUTE-ONCE-FILTER-MANY SWEEP
# ============================================================================

def compute_metrics_for_graph(G: nx.Graph, topk: int) -> dict:
    """Extract all metrics from a prebuilt graph."""
    n_total = G.number_of_nodes()
    edges = G.number_of_edges()
    nodes_connected = sum(1 for _, d in G.degree() if d > 0)

    density = nx.density(G) if n_total > 1 else float("nan")
    trans = nx.transitivity(G) if edges > 0 else float("nan")
    med_deg = float(np.median([d for _, d in G.degree()])) if n_total > 0 else float("nan")
    mpl = safe_mean_path_length(G)
    mod = modularity_louvain(G)
    cstats = component_stats(G)
    hub_list = topk_hubs(G, topk)

    return {
        "graph_nodes_total": n_total,
        "graph_nodes_connected": nodes_connected,
        "graph_edges": edges,
        "density": density,
        "transitivity": trans,
        "modularity_louvain": mod,
        **cstats,
        "median_degree": med_deg,
        "mean_path_length_LCC": mpl,
        "isolated_vip_features": cstats["n_singletons"],
        "_hub_list": hub_list,  # internal, stripped before CSV output
    }


def run_sweep(leaf_X_all, root_X_all, leaf_meta, root_meta,
              vip_leaf_df, vip_root_df,
              vip_thresholds, rho_thresholds, fdr_alpha, topk):
    """Compute-once-filter-many threshold sweep.

    For each (vip_thr, tissue, genotype) we compute the Spearman correlation
    and p-value matrices once. Then for each rho_thr we just filter the
    precomputed matrices -- this is a cheap mask operation. Eliminates
    len(rho_thresholds)-fold redundant Spearman recomputation and avoids
    multiprocessing pickle overhead on Windows.
    """
    all_metrics = []
    all_jacc = []
    total_combos = len(vip_thresholds) * len(rho_thresholds)
    done = 0

    for vip_thr in vip_thresholds:
        # Select VIP features for this threshold
        leaf_vips = set(vip_leaf_df.loc[vip_leaf_df["VIP_Score"] >= vip_thr, "Metabolite"].astype(str))
        root_vips = set(vip_root_df.loc[vip_root_df["VIP_Score"] >= vip_thr, "Metabolite"].astype(str))

        # Precompute Spearman matrices per (tissue, genotype) -- done ONCE per vip_thr
        precomputed = {}  # key: (tissue, genotype) -> (corr_df, p_df, used_cols, vip_total, n_all_na_dropped)

        for tissue, meta_df, X_all, vips in [
            ("Leaf", leaf_meta, leaf_X_all, leaf_vips),
            ("Root", root_meta, root_X_all, root_vips),
        ]:
            for genotype in ["G1", "G2"]:
                mask = (meta_df["Genotype"].astype(str) == genotype) & (meta_df["Treatment"] == 1)
                cols = [c for c in X_all.columns if c in vips]
                vip_total = len(cols)

                X = X_all.loc[mask, cols].copy()
                cols_before_na = X.shape[1]
                X = X.loc[:, X.notna().any(axis=0)]  # drop all-NA columns
                n_all_na_dropped = cols_before_na - X.shape[1]
                used_cols = list(X.columns)

                if X.shape[1] >= 2:
                    corr_df, p_df = spearman_corr_p(X)
                else:
                    corr_df, p_df = None, None

                precomputed[(tissue, genotype)] = (corr_df, p_df, used_cols, vip_total, n_all_na_dropped)

        print(f"  VIP >= {vip_thr}: Spearman matrices precomputed for 4 tissue x genotype slices")

        # Now sweep rho thresholds cheaply by filtering precomputed matrices
        for rho_thr in rho_thresholds:
            hubs_this_combo = {}

            for tissue, genotype in [("Leaf", "G1"), ("Leaf", "G2"),
                                     ("Root", "G1"), ("Root", "G2")]:
                corr_df, p_df, used_cols, vip_total, n_all_na_dropped = precomputed[(tissue, genotype)]

                if corr_df is not None:
                    G = build_graph_from_corr(corr_df, p_df, rho_thr=rho_thr, fdr_alpha=fdr_alpha)
                else:
                    G = nx.Graph()
                    G.add_nodes_from(used_cols)

                metrics = compute_metrics_for_graph(G, topk)
                hub_list = metrics.pop("_hub_list")
                hubs_this_combo[(tissue, genotype)] = hub_list

                row = {
                    "vip_threshold": vip_thr,
                    "rho_threshold": rho_thr,
                    "fdr_alpha": fdr_alpha,
                    "tissue": tissue,
                    "genotype": genotype,
                    "vip_selected_features": vip_total,
                    "used_features_after_na_drop": len(used_cols),
                    "n_all_na_dropped": n_all_na_dropped,
                    **metrics,
                }
                all_metrics.append(row)

            # Jaccard hub overlap per genotype
            for genotype in ["G1", "G2"]:
                leaf_h = hubs_this_combo.get(("Leaf", genotype), [])
                root_h = hubs_this_combo.get(("Root", genotype), [])
                all_jacc.append({
                    "vip_threshold": vip_thr,
                    "rho_threshold": rho_thr,
                    "genotype": genotype,
                    "topk": topk,
                    "hub_jaccard_leaf_vs_root": jaccard(leaf_h, root_h),
                    "leaf_hubs": ";".join(leaf_h),
                    "root_hubs": ";".join(root_h),
                })

            done += 1
            print(f"    [{done}/{total_combos}] VIP>={vip_thr}, |rho|>{rho_thr} -- done")

    return pd.DataFrame(all_metrics), pd.DataFrame(all_jacc)


# ============================================================================
# 3. SINGLETON-MISSINGNESS AUDIT (section 2.7, bugs fixed)
# ============================================================================

def singleton_missingness_audit(leaf_df, root_df, leaf_meta, root_meta,
                                vip_combined_set, rho_thr, fdr_alpha,
                                missing_s3b_df=None):
    """Cross-reference singleton/isolated features with missingness data.

    Fixes from v3 checklist:
      - Line 884: extra indent on miss = pd.read_csv(...)
      - Line 892: extra indent on for tissue_name, df in ...
      - Scalar guard for spearmanr when df has < 3 features (now handled
        in spearman_corr_p for >= 2 features; only skip if < 2)
      - Reports n_all_na_dropped to answer "are isolates just data-quality
        artifacts?"
    """
    rows = []

    for tissue_name, full_df, meta_df in [
        ("Leaf", leaf_df, leaf_meta),
        ("Root", root_df, root_meta),
    ]:
        X_all = full_df.iloc[:, 10:]  # metabolite columns

        # Filter to VIP features present in this tissue
        cols = [c for c in X_all.columns if c in vip_combined_set]

        for genotype in ["G1", "G2"]:
            mask = (meta_df["Genotype"].astype(str) == genotype) & (meta_df["Treatment"] == 1)
            X = X_all.loc[mask, cols].copy()
            cols_before = X.shape[1]
            X = X.dropna(axis=1, how="all")
            n_all_na_dropped = cols_before - X.shape[1]

            if X.shape[1] < 2:
                continue

            corr_df, p_df = spearman_corr_p(X)
            G = build_graph_from_corr(corr_df, p_df, rho_thr=rho_thr, fdr_alpha=fdr_alpha)
            cstats = component_stats(G)

            # Identify isolates
            isolates = [n for n, d in G.degree() if d == 0]
            connected = [n for n, d in G.degree() if d > 0]

            # Compute or look up missingness
            if missing_s3b_df is not None:
                miss_col = f"{tissue_name}_missing_pct"
                if miss_col in missing_s3b_df.columns:
                    miss = missing_s3b_df.set_index("Metabolite")
                    miss_iso = miss.loc[miss.index.isin(isolates), miss_col].dropna().values
                    miss_conn = miss.loc[miss.index.isin(connected), miss_col].dropna().values
                else:
                    per_feat_miss = (X.isna().mean(axis=0) * 100.0)
                    miss_iso = per_feat_miss.loc[per_feat_miss.index.isin(isolates)].values
                    miss_conn = per_feat_miss.loc[per_feat_miss.index.isin(connected)].values
            else:
                per_feat_miss = (X.isna().mean(axis=0) * 100.0)
                miss_iso = per_feat_miss.loc[per_feat_miss.index.isin(isolates)].values
                miss_conn = per_feat_miss.loc[per_feat_miss.index.isin(connected)].values

            row = {
                "tissue": tissue_name,
                "genotype": genotype,
                "rho_thr": rho_thr,
                "fdr_alpha": fdr_alpha,
                "n_vip_selected": cols_before,
                "n_all_na_dropped": n_all_na_dropped,
                "n_features_in_network": G.number_of_nodes(),
                "n_edges": G.number_of_edges(),
                **cstats,
                "n_isolates": len(isolates),
                "isolates_median_missing_pct": float(np.median(miss_iso)) if len(miss_iso) > 0 else float("nan"),
                "connected_median_missing_pct": float(np.median(miss_conn)) if len(miss_conn) > 0 else float("nan"),
                "isolates_mean_missing_pct": float(np.mean(miss_iso)) if len(miss_iso) > 0 else float("nan"),
                "connected_mean_missing_pct": float(np.mean(miss_conn)) if len(miss_conn) > 0 else float("nan"),
                "isolates_frac_missing_ge30": float(np.mean(miss_iso >= 30.0)) if len(miss_iso) > 0 else float("nan"),
                "connected_frac_missing_ge30": float(np.mean(miss_conn >= 30.0)) if len(miss_conn) > 0 else float("nan"),
                "isolate_names": ";".join(sorted(isolates)),
            }
            rows.append(row)

    return pd.DataFrame(rows)


# ============================================================================
# 4. FIGURE GENERATION (section 2.2)
# ============================================================================

def make_figures(metrics_df, jacc_df, out_dir, topk):
    """Generate Supp Figs S1 and S2."""
    if not HAS_MPL:
        print("[SKIP] matplotlib not available; skipping figures.")
        return

    # ── Style map: distinct markers and colors per tissue x genotype ──
    style_map = {
        ("Leaf", "G1"):  {"color": "#2ca02c", "marker": "o", "ls": "-"},
        ("Leaf", "G2"):  {"color": "#98df8a", "marker": "s", "ls": "--"},
        ("Root", "G1"):  {"color": "#d62728", "marker": "^", "ls": "-"},
        ("Root", "G2"):  {"color": "#ff9896", "marker": "D", "ls": "--"},
    }

    # ── Supp Fig S1: Threshold sensitivity at primary VIP ──
    primary_vip = 1.0
    s1_data = metrics_df[metrics_df["vip_threshold"] == primary_vip].copy()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f"Supplementary Figure S1: Network robustness across |rho| thresholds (VIP >= {primary_vip})",
                 fontsize=11, fontweight="bold", y=1.02)

    for (tissue, genotype), df in s1_data.groupby(["tissue", "genotype"]):
        sty = style_map[(tissue, genotype)]
        df = df.sort_values("rho_threshold")
        label = f"{tissue} {genotype}"

        axes[0].plot(df["rho_threshold"], df["density"],
                     marker=sty["marker"], color=sty["color"], ls=sty["ls"], label=label, linewidth=2)
        axes[1].plot(df["rho_threshold"], df["modularity_louvain"],
                     marker=sty["marker"], color=sty["color"], ls=sty["ls"], label=label, linewidth=2)

    axes[0].set_title("Network Density", fontsize=10)
    axes[0].set_xlabel("|rho| threshold")
    axes[0].set_ylabel("Density")
    axes[0].axvline(x=0.7, color="gray", ls=":", alpha=0.6, label="Primary threshold")

    axes[1].set_title("Louvain Modularity", fontsize=10)
    axes[1].set_xlabel("|rho| threshold")
    axes[1].set_ylabel("Modularity Q")
    axes[1].axvline(x=0.7, color="gray", ls=":", alpha=0.6)

    # Hub Jaccard
    jacc_primary = jacc_df[jacc_df["vip_threshold"] == primary_vip].copy()
    for genotype, df in jacc_primary.groupby("genotype"):
        df = df.sort_values("rho_threshold")
        color = "#1f77b4" if genotype == "G1" else "#aec7e8"
        marker = "o" if genotype == "G1" else "s"
        ls = "-" if genotype == "G1" else "--"
        axes[2].plot(df["rho_threshold"], df["hub_jaccard_leaf_vs_root"],
                     marker=marker, color=color, ls=ls, label=f"{genotype}", linewidth=2)

    axes[2].set_title(f"Leaf-Root Hub Overlap (top {topk})", fontsize=10)
    axes[2].set_xlabel("|rho| threshold")
    axes[2].set_ylabel("Jaccard index")
    axes[2].axvline(x=0.7, color="gray", ls=":", alpha=0.6)

    # Shared legend
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        for hi, li in zip(h, l):
            if li not in labels:
                handles.append(hi)
                labels.append(li)
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               bbox_to_anchor=(0.5, -0.08), frameon=True, fontsize=9)

    fig.tight_layout()
    out_s1 = out_dir / "Supp_Fig_S1_ThresholdSensitivity.pdf"
    fig.savefig(out_s1, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  [OK] {out_s1}")

    # ── Supp Fig S2: VIP sensitivity at primary rho=0.7 ──
    rho_default = 0.7
    s2_data = metrics_df[metrics_df["rho_threshold"] == rho_default].copy()

    if s2_data.empty:
        print(f"  [WARN] No data at rho={rho_default} for VIP sensitivity figure. Skipping S2.")
        return

    fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
    fig2.suptitle(f"Supplementary Figure S2: VIP threshold sensitivity at |rho| = {rho_default}",
                  fontsize=11, fontweight="bold", y=1.02)

    metric_panels = [
        ("density", "Density", axes2[0]),
        ("modularity_louvain", "Modularity Q", axes2[1]),
        ("transitivity", "Transitivity", axes2[2]),
    ]

    for metric_col, ylabel, ax in metric_panels:
        for (tissue, genotype), df in s2_data.groupby(["tissue", "genotype"]):
            sty = style_map[(tissue, genotype)]
            df = df.sort_values("vip_threshold")
            ax.plot(df["vip_threshold"], df[metric_col],
                    marker=sty["marker"], color=sty["color"], ls=sty["ls"],
                    label=f"{tissue} {genotype}", linewidth=2)
        ax.set_xlabel("VIP threshold")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel, fontsize=10)
        ax.axvline(x=1.0, color="gray", ls=":", alpha=0.6, label="Primary VIP")

    handles, labels = [], []
    for ax in axes2:
        h, l = ax.get_legend_handles_labels()
        for hi, li in zip(h, l):
            if li not in labels:
                handles.append(hi)
                labels.append(li)
    fig2.legend(handles, labels, loc="lower center", ncol=len(labels),
                bbox_to_anchor=(0.5, -0.08), frameon=True, fontsize=9)

    fig2.tight_layout()
    out_s2 = out_dir / "Supp_Fig_S2_VIPSensitivity.pdf"
    fig2.savefig(out_s2, bbox_inches="tight", dpi=300)
    plt.close(fig2)
    print(f"  [OK] {out_s2}")


# ============================================================================
# 5. VIP FILE LOADER
# ============================================================================

def load_vip_file(path: Path) -> pd.DataFrame:
    """Load VIP scores with flexible column name detection."""
    df = pd.read_csv(path)

    # Normalise first column to "Metabolite"
    df = df.rename(columns={df.columns[0]: "Metabolite"})

    # Find VIP score column
    if "VIP_Score" not in df.columns:
        for c in df.columns:
            if "vip" in c.lower() and c != "Metabolite":
                df = df.rename(columns={c: "VIP_Score"})
                break

    if "VIP_Score" not in df.columns:
        if len(df.columns) >= 2:
            df = df.rename(columns={df.columns[1]: "VIP_Score"})
        else:
            raise ValueError(f"Cannot find VIP score column in {path}. Columns: {list(df.columns)}")

    df["Metabolite"] = df["Metabolite"].astype(str)
    df["VIP_Score"] = pd.to_numeric(df["VIP_Score"], errors="coerce")
    return df[["Metabolite", "VIP_Score"]].dropna()


# ============================================================================
# 6. MAIN
# ============================================================================

def main():
    t0 = time.time()

    ap = argparse.ArgumentParser(
        description="Network robustness analysis for TPC-2025-1173 revision",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    ap.add_argument("--data-root", required=True, help="Project data root directory")
    ap.add_argument("--out", required=True, help="Output directory for all results")

    # Data file relative paths
    ap.add_argument("--leaf-data-rel", default=os.path.join("data", "data", "n_p_l.csv"))
    ap.add_argument("--root-data-rel", default=os.path.join("data", "data", "n_p_r.csv"))
    ap.add_argument("--vip-leaf-rel", default=os.path.join("output", "results", "updated_PLS", "vip_scores_all_leaf.csv"))
    ap.add_argument("--vip-root-rel", default=os.path.join("output", "results", "updated_PLS", "vip_scores_all_root.csv"))
    ap.add_argument("--vip-combined-rel",
                    default=os.path.join("output", "results", "spearman",
                                         "VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv"),
                    help="Combined VIP file for singleton audit (features with VIP > 1.0)")

    # Analysis parameters
    ap.add_argument("--vip-thresholds", default="1.0,0.8",
                    help="Comma-separated VIP thresholds to sweep")
    ap.add_argument("--rho-thresholds", default="0.6,0.65,0.7,0.75,0.8",
                    help="Comma-separated |rho| thresholds to sweep")
    ap.add_argument("--fdr-alpha", type=float, default=0.05)
    ap.add_argument("--topk", type=int, default=20, help="Top-k hubs for Jaccard overlap")

    # Options
    ap.add_argument("--no-plot", action="store_true", help="Skip figure generation")
    ap.add_argument("--missing-s3b", default=None,
                    help="Path to Supp_Table_S3b_MissingnessPerFeature.csv (optional)")

    args = ap.parse_args()

    data_root = Path(args.data_root)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    vip_thresholds = sorted([float(x.strip()) for x in args.vip_thresholds.split(",")], reverse=True)
    rho_thresholds = sorted([float(x.strip()) for x in args.rho_thresholds.split(",")])

    print("=" * 70)
    print("NETWORK ROBUSTNESS ANALYSIS v2 -- TPC-2025-1173 Revision")
    print("=" * 70)
    print(f"  Data root:       {data_root}")
    print(f"  Output dir:      {out_dir}")
    print(f"  VIP thresholds:  {vip_thresholds}")
    print(f"  |rho| thresholds: {rho_thresholds}")
    print(f"  FDR alpha:       {args.fdr_alpha}")
    print(f"  Top-k hubs:      {args.topk}")
    print(f"  Architecture:    compute-once-filter-many (no multiprocessing overhead)")
    print(f"  Figures:         {'Yes' if not args.no_plot else 'No'}")
    print(f"  Sample filter:   Treatment=1 (stressed only, matching manuscript)")
    print()

    # ── Load data ──
    print("[1/5] Loading data...")

    leaf_path = data_root / args.leaf_data_rel
    root_path = data_root / args.root_data_rel

    if not leaf_path.exists():
        sys.exit(f"[ERROR] Leaf data not found: {leaf_path}")
    if not root_path.exists():
        sys.exit(f"[ERROR] Root data not found: {root_path}")

    leaf = pd.read_csv(leaf_path)
    root = pd.read_csv(root_path)

    leaf_meta = leaf.iloc[:, :10]
    root_meta = root.iloc[:, :10]
    leaf_X_all = leaf.iloc[:, 10:]
    root_X_all = root.iloc[:, 10:]

    print(f"  Leaf: {leaf.shape[0]} samples x {leaf_X_all.shape[1]} features")
    print(f"  Root: {root.shape[0]} samples x {root_X_all.shape[1]} features")

    # ── Load VIP files with safety check for R2-4 ──
    vip_leaf_path = data_root / args.vip_leaf_rel
    vip_root_path = data_root / args.vip_root_rel
    has_per_tissue_vip = vip_leaf_path.exists() and vip_root_path.exists()

    min_vip_thr = min(vip_thresholds)

    if not has_per_tissue_vip and min_vip_thr < 1.0:
        # CRITICAL: fallback to combined VIP (already filtered to >1.0) would make
        # VIP=0.8 identical to VIP=1.0, silently invalidating the R2-4 analysis.
        sys.exit(
            f"[ERROR] Per-tissue VIP files not found:\n"
            f"  Leaf: {vip_leaf_path}\n"
            f"  Root: {vip_root_path}\n"
            f"\n"
            f"You requested VIP thresholds {vip_thresholds} which includes values < 1.0.\n"
            f"The combined VIP file is already filtered to VIP > 1.0, so using it as a\n"
            f"fallback would make VIP=0.8 and VIP=1.0 produce identical results -- silently\n"
            f"defeating the R2-4 sensitivity analysis.\n"
            f"\n"
            f"Fix: Re-run pls_tissue.py to output ALL VIP scores (unfiltered), then point\n"
            f"--vip-leaf-rel and --vip-root-rel to those files.\n"
            f"\n"
            f"  python .\\2_wheat_analysis\\pls_tissue.py\n"
            f"\n"
            f"See HANDOVER item #4: 'Does pls_tissue.py output all VIP scores or only >1.0?'"
        )

    if has_per_tissue_vip:
        vip_leaf_df = load_vip_file(vip_leaf_path)
        vip_root_df = load_vip_file(vip_root_path)
        print(f"  VIP Leaf: {len(vip_leaf_df)} features (full, unfiltered)")
        print(f"  VIP Root: {len(vip_root_df)} features (full, unfiltered)")

        # Sanity: verify there are actually features below 1.0
        if min_vip_thr < 1.0:
            n_below_leaf = (vip_leaf_df["VIP_Score"] < 1.0).sum()
            n_below_root = (vip_root_df["VIP_Score"] < 1.0).sum()
            if n_below_leaf == 0 and n_below_root == 0:
                print(f"  [WARN] VIP files contain NO features with VIP < 1.0.")
                print(f"         VIP sensitivity test may be trivial. Check pls_tissue.py output.")
            else:
                print(f"  Features with VIP < 1.0: Leaf={n_below_leaf}, Root={n_below_root}")
    else:
        # All thresholds are >= 1.0, safe to use combined file as fallback
        vip_combined_path = data_root / args.vip_combined_rel
        if not vip_combined_path.exists():
            sys.exit(f"[ERROR] Combined VIP file not found: {vip_combined_path}")
        vip_leaf_df = load_vip_file(vip_combined_path)
        vip_root_df = vip_leaf_df.copy()
        print(f"  VIP (combined fallback): {len(vip_leaf_df)} features")
        print(f"  [INFO] All VIP thresholds >= 1.0, combined fallback is safe.")

    # Load combined VIP for singleton audit
    vip_combined_path = data_root / args.vip_combined_rel
    if vip_combined_path.exists():
        vip_combined_df = pd.read_csv(vip_combined_path)
        vip_combined_df = vip_combined_df.rename(columns={vip_combined_df.columns[0]: "Metabolite"})
        vip_combined_set = set(vip_combined_df["Metabolite"].astype(str))
        print(f"  VIP combined (>1.0): {len(vip_combined_set)} features")
    else:
        print(f"  [WARN] Combined VIP file not found: {vip_combined_path}")
        vip_combined_set = set(vip_leaf_df.loc[vip_leaf_df["VIP_Score"] >= 1.0, "Metabolite"])
        print(f"  Using leaf VIP >= 1.0 as fallback: {len(vip_combined_set)} features")

    # Load missingness data if provided
    missing_s3b_df = None
    if args.missing_s3b and Path(args.missing_s3b).exists():
        missing_s3b_df = pd.read_csv(args.missing_s3b)
        if "Metabolite" not in missing_s3b_df.columns:
            missing_s3b_df = missing_s3b_df.rename(columns={missing_s3b_df.columns[0]: "Metabolite"})
        print(f"  Missingness file: {len(missing_s3b_df)} features loaded")

    print()

    # ── Threshold sweep ──
    n_spearman = len(vip_thresholds) * 4  # VIP x (2 tissues x 2 genotypes)
    n_filter = len(vip_thresholds) * len(rho_thresholds)
    print(f"[2/5] Running threshold sensitivity sweep...")
    print(f"  Spearman computations: {n_spearman} (once per VIP x tissue x genotype)")
    print(f"  Filtering passes:      {n_filter} (VIP x rho combos)")
    print(f"  Total networks built:  {n_filter * 4}")

    metrics_df, jacc_df = run_sweep(
        leaf_X_all=leaf_X_all, root_X_all=root_X_all,
        leaf_meta=leaf_meta, root_meta=root_meta,
        vip_leaf_df=vip_leaf_df, vip_root_df=vip_root_df,
        vip_thresholds=vip_thresholds,
        rho_thresholds=rho_thresholds,
        fdr_alpha=args.fdr_alpha,
        topk=args.topk,
    )

    # Save Supp Table S1
    out_s1 = out_dir / "Supp_Table_S1_WheatNetworkSensitivity_Metrics.csv"
    metrics_df.to_csv(out_s1, index=False)
    print(f"  [OK] {out_s1}")
    print(f"       {len(metrics_df)} rows: {len(vip_thresholds)} VIP x {len(rho_thresholds)} rho x 2 tissues x 2 genotypes")

    # Save Supp Table S2
    out_s2 = out_dir / "Supp_Table_S2_WheatNetworkSensitivity_HubJaccard.csv"
    jacc_df.to_csv(out_s2, index=False)
    print(f"  [OK] {out_s2}")
    print()

    # ── Sanity check: does the asymmetry hold? ──
    print("[3/5] Sanity check: leaf-root asymmetry across thresholds...")
    primary_rho = 0.7
    primary_vip_check = 1.0
    primary = metrics_df[
        (metrics_df["vip_threshold"] == primary_vip_check) &
        (metrics_df["rho_threshold"] == primary_rho)
    ]

    if not primary.empty:
        for genotype in ["G1", "G2"]:
            leaf_row = primary[(primary["tissue"] == "Leaf") & (primary["genotype"] == genotype)]
            root_row = primary[(primary["tissue"] == "Root") & (primary["genotype"] == genotype)]
            if not leaf_row.empty and not root_row.empty:
                ld = leaf_row.iloc[0]["density"]
                rd = root_row.iloc[0]["density"]
                lm = leaf_row.iloc[0]["modularity_louvain"]
                rm = root_row.iloc[0]["modularity_louvain"]
                asymmetry_ok = (ld > rd) and (rm > lm)
                status = "CONFIRMED" if asymmetry_ok else "!! CHECK !!"
                print(f"  {genotype}: Leaf density={ld:.4f} > Root density={rd:.4f}, "
                      f"Root mod={rm:.3f} > Leaf mod={lm:.3f} -> {status}")

    # Check across all thresholds
    asymmetry_density_count = 0
    asymmetry_mod_count = 0
    total_checks = 0
    for _, row in metrics_df.iterrows():
        if row["tissue"] == "Leaf":
            partner = metrics_df[
                (metrics_df["vip_threshold"] == row["vip_threshold"]) &
                (metrics_df["rho_threshold"] == row["rho_threshold"]) &
                (metrics_df["tissue"] == "Root") &
                (metrics_df["genotype"] == row["genotype"])
            ]
            if not partner.empty:
                total_checks += 1
                if row["density"] > partner.iloc[0]["density"]:
                    asymmetry_density_count += 1
                if row["modularity_louvain"] < partner.iloc[0]["modularity_louvain"]:
                    asymmetry_mod_count += 1

    if total_checks > 0:
        pct_d = (asymmetry_density_count / total_checks) * 100
        pct_m = (asymmetry_mod_count / total_checks) * 100
        print(f"  Leaf > Root density:      {asymmetry_density_count}/{total_checks} ({pct_d:.0f}%)")
        print(f"  Root > Leaf modularity:   {asymmetry_mod_count}/{total_checks} ({pct_m:.0f}%)")
    print()

    # ── Singleton-missingness audit ──
    print("[4/5] Running singleton-missingness audit (R1-Meth-6)...")
    audit_df = singleton_missingness_audit(
        leaf_df=leaf, root_df=root,
        leaf_meta=leaf_meta, root_meta=root_meta,
        vip_combined_set=vip_combined_set,
        rho_thr=0.7, fdr_alpha=args.fdr_alpha,
        missing_s3b_df=missing_s3b_df,
    )

    out_s7 = out_dir / "Supp_Table_S7_WheatComponentSingletonAudit.csv"
    audit_df.to_csv(out_s7, index=False)
    print(f"  [OK] {out_s7}")

    if not audit_df.empty:
        for _, row in audit_df.iterrows():
            iso_miss = row.get("isolates_median_missing_pct", float("nan"))
            conn_miss = row.get("connected_median_missing_pct", float("nan"))
            na_drop = row.get("n_all_na_dropped", 0)
            print(f"  {row['tissue']} {row['genotype']}: "
                  f"{row['n_isolates']} isolates (median miss={iso_miss:.1f}%) vs "
                  f"connected (median miss={conn_miss:.1f}%), "
                  f"LCC covers {row['lcc_frac']:.1%} of nodes, "
                  f"{na_drop} features dropped (all-NA)")
    print()

    # ── Figures ──
    if not args.no_plot:
        print("[5/5] Generating supplementary figures...")
        make_figures(metrics_df, jacc_df, out_dir, topk=args.topk)
    else:
        print("[5/5] Figure generation skipped (--no-plot).")

    elapsed = time.time() - t0
    print()
    print("=" * 70)
    print(f"DONE in {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"All outputs in: {out_dir}")
    print()
    print("Output files:")
    print(f"  Supp Table S1: {out_s1.name}")
    print(f"  Supp Table S2: {out_s2.name}")
    print(f"  Supp Table S7: {out_s7.name}")
    if not args.no_plot:
        print(f"  Supp Fig  S1:  Supp_Fig_S1_ThresholdSensitivity.pdf")
        print(f"  Supp Fig  S2:  Supp_Fig_S2_VIPSensitivity.pdf")
    print("=" * 70)


if __name__ == "__main__":
    main()