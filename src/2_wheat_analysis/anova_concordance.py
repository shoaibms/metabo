#!/usr/bin/env python3
"""
anova_concordance.py — ANOVA concordance check for TPC-2025-1173 revision.

Addresses:
  R1-Meth-4  "Why Cliff's δ not ANOVA?"

Purpose:
  Defensive supplement showing that for the subset of features satisfying
  ANOVA's normality assumption, parametric and non-parametric tests yield
  concordant significance calls.  This confirms the choice of Cliff's δ
  does not bias biological conclusions.

Approach:
  For each tissue × genotype × day combination:
    1. Separate stressed (Treatment=1) vs control (Treatment=0) groups.
    2. Per metabolite, run Shapiro-Wilk in both groups.
    3. Features passing normality in BOTH groups (P > 0.05) enter the
       concordance test.
    4. For those features: run one-way ANOVA (F-test) AND Mann-Whitney U.
    5. Apply BH-FDR correction within each comparison block.
    6. Compare significance calls (concordant / discordant at FDR < 0.05).

Note on Treatment filter:
  Other revision scripts filter to Treatment=1 only (stressed) for network
  construction. This script necessarily uses BOTH treatment groups because
  the comparison IS the treatment effect (stressed vs control), matching the
  original stat_tests.py (FILE 31) analysis that the reviewer questioned.

Outputs:
  Supp_Table_SD_ANOVAConcordance_Detail.csv  — per-feature results (GitHub only; not in journal submission)
  Supp_Table_SD_ANOVAConcordance_Summary.csv — aggregate concordance stats (submission: Supp Table S9)

Usage:
  python anova_concordance.py --data-root "<data_root>" --out "output/revision_pack"

  Optional:
    --leaf-data-rel  (default: data\\data\\n_p_l.csv)
    --root-data-rel  (default: data\\data\\n_p_r.csv)
    --fdr-alpha      (default: 0.05)
    --shapiro-alpha  (default: 0.05)
"""

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

METADATA_COLS = 10  # First 10 columns are metadata in n_p_l/r.csv


# ============================================================================
# BH-FDR (same implementation as other revision scripts)
# ============================================================================

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR q-values (NaN-safe)."""
    p = np.asarray(pvals, dtype=float)
    p = np.nan_to_num(p, nan=1.0, posinf=1.0, neginf=1.0)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out


# ============================================================================
# CORE ANALYSIS
# ============================================================================

def run_concordance(df: pd.DataFrame, tissue: str,
                    shapiro_alpha: float, fdr_alpha: float) -> pd.DataFrame:
    """Run ANOVA vs Mann-Whitney concordance for one tissue.

    For each genotype × day combination, compares Treatment=1 vs Treatment=0
    per metabolite. Only features passing Shapiro-Wilk in both groups are tested.

    Parameters
    ----------
    df : pd.DataFrame
        Full data matrix (metadata + features).
    tissue : str
        Label ("Leaf" or "Root").
    shapiro_alpha : float
        Significance threshold for Shapiro-Wilk normality test.
    fdr_alpha : float
        FDR threshold for significance calls.

    Returns
    -------
    pd.DataFrame with per-feature results.
    """
    X = df.iloc[:, METADATA_COLS:]
    meta = df.iloc[:, :METADATA_COLS]

    # Identify metabolite columns
    feat_cols = list(X.columns)

    genotypes = sorted(meta["Genotype"].astype(str).unique())
    days = sorted(meta["Day"].unique())

    rows = []

    for genotype in genotypes:
        for day in days:
            # Get stressed and control groups
            mask_s = ((meta["Genotype"].astype(str) == genotype) &
                      (meta["Treatment"] == 1) &
                      (meta["Day"] == day))
            mask_c = ((meta["Genotype"].astype(str) == genotype) &
                      (meta["Treatment"] == 0) &
                      (meta["Day"] == day))

            stressed = X.loc[mask_s]
            control = X.loc[mask_c]

            n_s = stressed.shape[0]
            n_c = control.shape[0]

            if n_s < 3 or n_c < 3:
                continue

            # Per-feature tests
            block_rows = []
            for feat in feat_cols:
                s_vals = stressed[feat].dropna().values
                c_vals = control[feat].dropna().values

                if len(s_vals) < 3 or len(c_vals) < 3:
                    continue

                # Shapiro-Wilk normality
                try:
                    _, sw_p_s = sp_stats.shapiro(s_vals)
                    _, sw_p_c = sp_stats.shapiro(c_vals)
                except Exception:
                    continue

                passes_normality = (sw_p_s > shapiro_alpha) and (sw_p_c > shapiro_alpha)

                # Mann-Whitney U (always computed)
                try:
                    mw_stat, mw_p = sp_stats.mannwhitneyu(
                        s_vals, c_vals, alternative="two-sided"
                    )
                except Exception:
                    mw_stat, mw_p = np.nan, np.nan

                # One-way ANOVA (always computed for normal features)
                if passes_normality:
                    try:
                        f_stat, anova_p = sp_stats.f_oneway(s_vals, c_vals)
                    except Exception:
                        f_stat, anova_p = np.nan, np.nan
                else:
                    f_stat, anova_p = np.nan, np.nan

                block_rows.append({
                    "tissue": tissue,
                    "genotype": genotype,
                    "day": day,
                    "feature": feat,
                    "n_stressed": len(s_vals),
                    "n_control": len(c_vals),
                    "shapiro_p_stressed": round(sw_p_s, 6),
                    "shapiro_p_control": round(sw_p_c, 6),
                    "passes_normality": passes_normality,
                    "mw_p_raw": mw_p,
                    "anova_p_raw": anova_p,
                    "anova_F": f_stat,
                })

            if not block_rows:
                continue

            # BH-FDR within this genotype × day block
            block_df = pd.DataFrame(block_rows)

            # FDR for Mann-Whitney (all features — for general reporting)
            mw_pvals = block_df["mw_p_raw"].values
            block_df["mw_q"] = bh_fdr(mw_pvals)

            # FDR for ANOVA (only normal features; others stay NaN)
            normal_mask = block_df["passes_normality"].values
            anova_q = np.full(len(block_df), np.nan)
            if normal_mask.sum() > 0:
                anova_raw = block_df.loc[normal_mask, "anova_p_raw"].values
                anova_q[normal_mask] = bh_fdr(anova_raw)
            block_df["anova_q"] = anova_q

            # Fair MW FDR on same normal-only subset (for concordance only)
            # This ensures both tests are corrected over identical test counts,
            # avoiding mechanical inflation of discordance from unequal FDR pools.
            mw_q_normal = np.full(len(block_df), np.nan)
            if normal_mask.sum() > 0:
                mw_raw_normal = block_df.loc[normal_mask, "mw_p_raw"].values
                mw_q_normal[normal_mask] = bh_fdr(mw_raw_normal)
            block_df["mw_q_normal"] = mw_q_normal

            # Significance calls
            block_df["mw_sig"] = block_df["mw_q"] < fdr_alpha
            block_df["anova_sig"] = block_df["anova_q"] < fdr_alpha
            # Fair comparison: MW significance on same subset as ANOVA
            block_df["mw_sig_normal"] = block_df["mw_q_normal"] < fdr_alpha

            # Concordance (only for features that pass normality)
            # Uses mw_sig_normal for fair apple-to-apple comparison
            block_df["concordant"] = np.where(
                block_df["passes_normality"],
                block_df["mw_sig_normal"] == block_df["anova_sig"],
                np.nan
            )

            rows.append(block_df)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, ignore_index=True)


def summarise_concordance(detail_df: pd.DataFrame) -> pd.DataFrame:
    """Compute aggregate concordance statistics.

    Returns one row per tissue (+ overall) with:
      - total features tested
      - n and % passing normality
      - n concordant / discordant among normal features
      - concordance percentage
    """
    if detail_df.empty:
        return pd.DataFrame()

    summaries = []
    groups = list(detail_df.groupby("tissue")) + [("Overall", detail_df)]

    for label, sub in groups:
        n_total = len(sub)
        n_normal = int(sub["passes_normality"].sum())
        pct_normal = round(n_normal / n_total * 100, 2) if n_total else 0

        normal_sub = sub[sub["passes_normality"]]
        n_concordant = int(normal_sub["concordant"].sum()) if len(normal_sub) else 0
        n_discordant = len(normal_sub) - n_concordant
        pct_concordant = round(n_concordant / len(normal_sub) * 100, 2) if len(normal_sub) else 0

        # Break down discordant cases
        if len(normal_sub) > 0:
            discord = normal_sub[~normal_sub["concordant"].astype(bool)]
            n_anova_only = int((discord["anova_sig"] & ~discord["mw_sig_normal"]).sum())
            n_mw_only = int((discord["mw_sig_normal"] & ~discord["anova_sig"]).sum())
        else:
            n_anova_only = 0
            n_mw_only = 0

        summaries.append({
            "scope": label,
            "n_comparisons_total": n_total,
            "n_passing_normality": n_normal,
            "pct_passing_normality": pct_normal,
            "n_concordant": n_concordant,
            "n_discordant": n_discordant,
            "n_discordant_anova_only_sig": n_anova_only,
            "n_discordant_mw_only_sig": n_mw_only,
            "pct_concordant": pct_concordant,
        })

    return pd.DataFrame(summaries)


# ============================================================================
# MAIN
# ============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="ANOVA vs Mann-Whitney concordance check (R1-Meth-4)."
    )
    ap.add_argument("--data-root", required=True,
                    help="Root data directory")
    ap.add_argument("--out", required=True,
                    help="Output directory")
    ap.add_argument("--leaf-data-rel", default=r"data\data\n_p_l.csv",
                    help="Relative path to leaf data CSV")
    ap.add_argument("--root-data-rel", default=r"data\data\n_p_r.csv",
                    help="Relative path to root data CSV")
    ap.add_argument("--fdr-alpha", type=float, default=0.05,
                    help="FDR significance threshold (default: 0.05)")
    ap.add_argument("--shapiro-alpha", type=float, default=0.05,
                    help="Shapiro-Wilk normality threshold (default: 0.05)")
    args = ap.parse_args()

    data_root = Path(args.data_root)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Normalise backslash defaults for cross-platform compatibility
    args.leaf_data_rel = args.leaf_data_rel.replace("\\", "/")
    args.root_data_rel = args.root_data_rel.replace("\\", "/")

    if not data_root.exists():
        sys.exit(f"[ERROR] Data root not found: {data_root}")

    t0 = time.time()

    print("=" * 70)
    print("anova_concordance.py — TPC-2025-1173 Revision (R1-Meth-4)")
    print("=" * 70)
    print(f"  Data root:      {data_root}")
    print(f"  Output dir:     {out_dir}")
    print(f"  FDR alpha:      {args.fdr_alpha}")
    print(f"  Shapiro alpha:  {args.shapiro_alpha}")
    print()

    # ── Load data ──
    leaf_path = data_root / args.leaf_data_rel
    root_path = data_root / args.root_data_rel

    if not leaf_path.exists():
        sys.exit(f"[ERROR] Leaf data not found: {leaf_path}")
    if not root_path.exists():
        sys.exit(f"[ERROR] Root data not found: {root_path}")

    leaf = pd.read_csv(leaf_path)
    root = pd.read_csv(root_path)

    print(f"  Leaf: {leaf.shape[0]} rows × {leaf.shape[1] - METADATA_COLS} features")
    print(f"  Root: {root.shape[0]} rows × {root.shape[1] - METADATA_COLS} features")

    # Verify required columns
    for name, df in [("Leaf", leaf), ("Root", root)]:
        for col in ["Genotype", "Treatment", "Day"]:
            if col not in df.columns:
                sys.exit(f"[ERROR] '{col}' column not found in {name} data.")

    print()

    # ── Run concordance analysis ──
    print("[1/3] Running concordance analysis...")
    results = []
    for tissue_name, df in [("Leaf", leaf), ("Root", root)]:
        print(f"  {tissue_name}...", end=" ", flush=True)
        tissue_df = run_concordance(df, tissue_name, args.shapiro_alpha, args.fdr_alpha)
        if not tissue_df.empty:
            results.append(tissue_df)
            n_normal = tissue_df["passes_normality"].sum()
            n_total = len(tissue_df)
            print(f"{n_total:,} comparisons, {n_normal:,} pass normality "
                  f"({n_normal/n_total*100:.1f}%)")
        else:
            print("no valid comparisons found.")

    if not results:
        print("\n[ERROR] No results produced. Check data structure.")
        sys.exit(1)

    detail_df = pd.concat(results, ignore_index=True)

    # ── Save detail table ──
    print()
    print("[2/3] Saving results...")
    out_detail = out_dir / "Supp_Table_SD_ANOVAConcordance_Detail.csv"
    detail_df.to_csv(out_detail, index=False)
    print(f"  [OK] {out_detail}")
    print(f"       {len(detail_df):,} rows total")

    # ── Summary ──
    summary_df = summarise_concordance(detail_df)
    out_summary = out_dir / "Supp_Table_SD_ANOVAConcordance_Summary.csv"
    summary_df.to_csv(out_summary, index=False)
    print(f"  [OK] {out_summary}")

    # ── Print manuscript-ready results ──
    print()
    print("[3/3] Results summary:")
    print("-" * 60)
    for _, row in summary_df.iterrows():
        print(f"  {row['scope']}:")
        print(f"    Total comparisons:    {row['n_comparisons_total']:,}")
        print(f"    Passing normality:    {row['n_passing_normality']:,} "
              f"({row['pct_passing_normality']:.1f}%)")
        print(f"    Concordant calls:     {row['n_concordant']:,} / "
              f"{row['n_passing_normality']:,} "
              f"({row['pct_concordant']:.1f}%)")
        if row['n_discordant'] > 0:
            print(f"    Discordant:           {row['n_discordant']:,} "
                  f"(ANOVA-only sig: {row['n_discordant_anova_only_sig']}, "
                  f"MW-only sig: {row['n_discordant_mw_only_sig']})")
        print()

    # Manuscript one-liner
    overall = summary_df[summary_df["scope"] == "Overall"]
    if not overall.empty:
        o = overall.iloc[0]
        print("  MANUSCRIPT SENTENCE:")
        print(f"  \"For the {o['pct_passing_normality']:.0f}% of feature-level comparisons")
        print(f"   satisfying normality assumptions (Shapiro-Wilk P > {args.shapiro_alpha}),")
        print(f"   ANOVA yielded concordant significance calls in")
        print(f"   {o['pct_concordant']:.1f}% of comparisons (Supplementary Table SD).\"")
        print()

    elapsed = time.time() - t0
    print("=" * 70)
    print(f"DONE in {elapsed:.1f}s")
    print(f"All outputs in: {out_dir}")
    print()
    print("Output files:")
    print(f"  Supp Table SD (detail):  {out_detail.name}")
    print(f"  Supp Table SD (summary): {out_summary.name}")
    print("=" * 70)


if __name__ == "__main__":
    main()