# =============================================================================
#  figure_s9_statistical_robustness.py
#  Supplementary Figure S9 — Statistical Power and Network Stability
#  TPC-2025-1173
#
#  Script location : <project_root>/rivision_pack_1/
#  Output location : <data_root>/output/revision_pack/figures/
#
#  DATA SOURCES
#  ---------------------------------------------------------------
#  Panels A–B  (power histograms):
#    output/results/initial_stat/initial_stat6e/tissue_comparison.csv
#    output/results/initial_stat/initial_stat6e/genotype_comparison.csv
#      column used: "Power"
#      NOTE: temporal_analysis.csv does not exist — Panel C omitted.
#
#  Panels D–F  (bootstrap boxplots):
#    output/revision_pack/Supp_Table_S6_WheatNetworkMetricStability.csv
#      columns: tissue, genotype, resample, iter,
#               density_all, modularity_all, transitivity_lcc
#    output/revision_pack/Supp_Table_S6b_WheatBootstrapSummary.csv
#      columns: tissue, genotype, metric, baseline,
#               bootstrap_ci_lo, bootstrap_ci_hi, jackknife_cv
#
#  Panel E G2 modularity Δ CI recomputed from raw S6 bootstrap iterations.
# =============================================================================

import os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
warnings.filterwarnings("ignore")
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "axes.spines.top": True, "axes.spines.right": True,
})

# =============================================================================
#  PATHS — USER CONFIGURATION: update DATA_ROOT for your local environment
#  All other paths below are derived from DATA_ROOT automatically.
# =============================================================================

DATA_ROOT  = r"C:\Users\ms\Desktop\data_chem_3_10"
STAT_DIR   = os.path.join(DATA_ROOT,
                          r"output\results\initial_stat\initial_stat6e")
PACK_DIR   = os.path.join(DATA_ROOT, r"output\revision_pack")
OUTPUT_DIR = os.path.join(DATA_ROOT, r"output\revision_pack\figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
#  LOAD DATA — power histograms (Panels A–B)
# =============================================================================

# Temporal file exists in one of three initial_stat directories.
# Search all three and use the first one found.
def _find_temporal():
    stat_variants = [
        os.path.join(DATA_ROOT, r"output\results\initial_stat\initial_stat6e"),
        os.path.join(DATA_ROOT, r"output\results\initial_stat\initial_stat6e_II"),
        os.path.join(DATA_ROOT, r"output\results\initial_stat\initial_stat6e_III"),
    ]
    for d in stat_variants:
        p = os.path.join(d, "temporal_analysis.csv")
        if os.path.isfile(p):
            print(f"[OK] temporal_analysis.csv found in: {d}")
            return p
    return None

TEMPORAL_PATH = _find_temporal()

POWER_FILES = {
    "tissue":   os.path.join(STAT_DIR, "tissue_comparison.csv"),
    "genotype": os.path.join(STAT_DIR, "genotype_comparison.csv"),
}
if TEMPORAL_PATH:
    POWER_FILES["temporal"] = TEMPORAL_PATH
else:
    print("[WARN] temporal_analysis.csv not found in any of the three expected locations.")

power_series = {}
for label, path in POWER_FILES.items():
    if not os.path.isfile(path):
        sys.exit(
            f"[ERROR] Power file not found:\n  {path}\n"
            f"Run export_revision_stats.py first (or check initial_stat directory name)."
        )
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    if "Power" not in df.columns:
        sys.exit(
            f"[ERROR] 'Power' column missing in {os.path.basename(path)}.\n"
            f"  Available columns: {list(df.columns)}"
        )
    pwr = df["Power"].dropna().astype(float)
    power_series[label] = pwr
    pct_80 = (pwr >= 0.80).mean() * 100
    print(f"[OK] {label}: {len(pwr):,} tests  |  median={pwr.median():.3f}  "
          f"|  >=0.80: {pct_80:.1f}%")

# =============================================================================
#  LOAD DATA — bootstrap stability (Panels D–F)
# =============================================================================

raw_path = os.path.join(PACK_DIR, "Supp_Table_S6_WheatNetworkMetricStability.csv")
sum_path = os.path.join(PACK_DIR, "Supp_Table_S6b_WheatBootstrapSummary.csv")

for p in [raw_path, sum_path]:
    if not os.path.isfile(p):
        sys.exit(
            f"[ERROR] Bootstrap file not found:\n  {p}\n"
            f"Run bootstrap_stability.py first."
        )

raw = pd.read_csv(raw_path)
raw.columns = [c.strip() for c in raw.columns]
summ = pd.read_csv(sum_path)
summ.columns = [c.strip() for c in summ.columns]

# Verify required columns
for col in ["tissue", "genotype", "resample", "iter",
            "density_all", "modularity_all", "transitivity_lcc"]:
    if col not in raw.columns:
        sys.exit(f"[ERROR] Expected column '{col}' not found in raw bootstrap file.\n"
                 f"  Available: {list(raw.columns)}")

print(f"[OK] Bootstrap raw: {len(raw)} rows  "
      f"({(raw['resample']=='bootstrap').sum()} bootstrap iterations)")
print(f"[OK] Bootstrap summary: {len(summ)} rows")

# Isolate bootstrap-only rows (exclude baseline and jackknife from boxplot data)
boot = raw[raw["resample"] == "bootstrap"].copy()

# Baseline values (all 7 replicates) — used as star marker
baseline = raw[raw["resample"] == "baseline"].copy()

# Jackknife CV — from summary file
def get_cv(tissue, genotype, metric):
    row = summ[(summ["tissue"] == tissue) &
               (summ["genotype"] == genotype) &
               (summ["metric"] == metric)]
    if row.empty:
        return float("nan")
    return float(row["jackknife_cv"].iloc[0])

def get_baseline(tissue, genotype, metric):
    row = baseline[(baseline["tissue"] == tissue) &
                   (baseline["genotype"] == genotype)]
    if row.empty or metric not in row.columns:
        return float("nan")
    return float(row[metric].iloc[0])

# =============================================================================
#  PANEL E — Compute G2 modularity leaf−root Δ from raw bootstrap iterations
#  (paired differences: same seed per genotype → same replicates per iter)
# =============================================================================

def compute_paired_delta(genotype, metric):
    """
    Compute paired bootstrap difference: Leaf − Root for given genotype+metric.
    Returns dict: baseline_diff, mean_diff, ci_lo, ci_hi, ci_excludes_zero
    """
    leaf_b = (boot[(boot["genotype"] == genotype) & (boot["tissue"] == "Leaf")]
              .sort_values("iter").reset_index(drop=True))
    root_b = (boot[(boot["genotype"] == genotype) & (boot["tissue"] == "Root")]
              .sort_values("iter").reset_index(drop=True))

    if leaf_b.empty or root_b.empty:
        return None

    n = min(len(leaf_b), len(root_b))
    diff = leaf_b[metric].values[:n] - root_b[metric].values[:n]
    diff = diff[~np.isnan(diff)]

    if len(diff) == 0:
        return None

    ci_lo = float(np.percentile(diff, 2.5))
    ci_hi = float(np.percentile(diff, 97.5))
    base_leaf = get_baseline(genotype, "Leaf", metric) if False else \
                get_baseline("Leaf", genotype, metric)
    base_root = get_baseline("Root", genotype, metric)
    base_diff = base_leaf - base_root if not (np.isnan(base_leaf) or np.isnan(base_root)) else float("nan")

    # For modularity: Leaf − Root expected to be negative (root > leaf)
    significant = ci_hi < 0

    return {
        "baseline_diff": base_diff,
        "mean_diff": float(np.mean(diff)),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "n_pairs": len(diff),
        "ci_excludes_zero": significant,
    }

g1_mod_delta = compute_paired_delta("G1", "modularity_all")
g2_mod_delta = compute_paired_delta("G2", "modularity_all")

if g2_mod_delta:
    sig_str = "significant" if g2_mod_delta["ci_excludes_zero"] else "not significant"
    print(f"[OK] G2 modularity Leaf-Root delta: mean={g2_mod_delta['mean_diff']:.3f}  "
          f"95% CI [{g2_mod_delta['ci_lo']:.3f}, {g2_mod_delta['ci_hi']:.3f}]  ({sig_str})")
if g1_mod_delta:
    print(f"[OK] G1 modularity Leaf-Root delta: mean={g1_mod_delta['mean_diff']:.3f}  "
          f"95% CI [{g1_mod_delta['ci_lo']:.3f}, {g1_mod_delta['ci_hi']:.3f}]")

print("\n[OK] All data loaded — building figure...\n")

# =============================================================================
#  PALETTE  (manuscript identity — from design brief + colours.py)
# =============================================================================

# Tissue × genotype (bootstrap boxes + histogram A)
STYLE = {
    ("Leaf", "G1"): "#3CB371",   # MediumSeaGreen
    ("Leaf", "G2"): "#2AA9B0",   # G2 leaf teal
    ("Root", "G1"): "#2E8B57",   # SeaGreen (darker)
    ("Root", "G2"): "#1B6F75",   # Darker teal
}
POWER_A_COL  = "#3CB371"   # tissue = primary comparison
POWER_B_COL  = "#BDBDBD"   # genotype/temporal = secondary, explicitly grey
CORAL        = "#DC6B4A"   # 0.80 threshold line + annotations
ANNO_BG      = "#FFFFFF"
ANNO_EC      = "#AAAAAA"
GRID_COL     = "#DDDDDD"
REF_COL      = "#999999"

TAG = 16; SUB = 14; AXT = 13; TCK = 11; LEG = 11; ANN = 10

# Group order for boxplots
GROUPS = [("Leaf", "G1"), ("Leaf", "G2"), ("Root", "G1"), ("Root", "G2")]
GROUP_LABELS = ["Leaf\nG1", "Leaf\nG2", "Root\nG1", "Root\nG2"]

# =============================================================================
#  HELPERS
# =============================================================================

def spine_style(ax, bg="white"):
    for s in ax.spines.values():
        s.set_linewidth(0.8); s.set_color("black")
    ax.yaxis.grid(True, alpha=0.28, lw=0.5, color=GRID_COL, ls="-")
    ax.xaxis.grid(False); ax.set_axisbelow(True)
    ax.tick_params(labelsize=TCK, length=3, width=0.6)
    ax.set_facecolor(bg)

def ptag(ax, tag, x=-0.14, y=1.06):
    ax.text(x, y, tag, transform=ax.transAxes,
            fontsize=TAG, fontweight="normal", va="bottom", ha="left",
            path_effects=[pe.withStroke(linewidth=1.2, foreground="white")])

def abox(ax, x, y, txt, ha="right", va="top", fs=ANN,
         bold=False, col="#333333", ec=ANNO_EC):
    ax.annotate(txt, xy=(x, y), xycoords=ax.transAxes,
                fontsize=fs, ha=ha, va=va, color=col,
                fontweight="bold" if bold else "normal",
                bbox=dict(boxstyle="round,pad=0.28", fc=ANNO_BG,
                          ec=ec, lw=0.7, alpha=0.95))

# =============================================================================
#  FIGURE LAYOUT  (asymmetric)
#
#  Top row (taller, ratio 1.25):
#    A — Tissue power (60%)  |  B — Genotype power (40%)
#
#  Bottom row (ratio 1.0):
#    D — Density (25%)  |  E — Modularity / hero (45%)  |  F — Transitivity (30%)
#
#  Total: 16" × 10"
# =============================================================================

fig = plt.figure(figsize=(13, 9), dpi=300, facecolor="white")
outer = gridspec.GridSpec(
    2, 1, figure=fig,
    hspace=0.74,
    height_ratios=[1.25, 1.0],
    top=0.92, bottom=0.14, left=0.07, right=0.97
)

# 3 top panels if temporal exists, else 2
_has_temporal = bool(TEMPORAL_PATH)
_top_ratios   = [0.42, 0.30, 0.28] if _has_temporal else [0.60, 0.40]
_top_ncols    = 3 if _has_temporal else 2
top = gridspec.GridSpecFromSubplotSpec(
    1, _top_ncols, subplot_spec=outer[0],
    wspace=0.32, width_ratios=_top_ratios
)
bot = gridspec.GridSpecFromSubplotSpec(
    1, 3, subplot_spec=outer[1],
    wspace=0.36, width_ratios=[0.25, 0.45, 0.30]
)

ax_A = fig.add_subplot(top[0])
ax_B = fig.add_subplot(top[1])
ax_C = fig.add_subplot(top[2]) if _has_temporal else None
ax_D = fig.add_subplot(bot[0])
ax_E = fig.add_subplot(bot[1])
ax_F = fig.add_subplot(bot[2])

# =============================================================================
#  POWER HISTOGRAM HELPER
# =============================================================================

def power_histogram(ax, pwr_series, label, fill_col, tag, subtitle=None, tag_x=-0.12):
    """
    Draws a single power distribution histogram.
    pwr_series : pd.Series of per-feature power values
    label      : comparison type string
    fill_col   : bar fill colour
    tag        : panel letter
    subtitle   : optional italic subtitle below panel title
    """
    spine_style(ax, bg="#FDFDFB")

    counts, edges = np.histogram(pwr_series.values, bins=30, range=(0, 1))
    bin_w  = edges[1] - edges[0]
    bin_cx = (edges[:-1] + edges[1:]) / 2

    # Gradient shade: lighter left, darker right
    cmax = counts.max() if counts.max() > 0 else 1
    for cx, cnt, lft, rgt in zip(bin_cx, counts, edges[:-1], edges[1:]):
        # bins right of threshold get slightly darker shade
        shade = fill_col
        alpha = 0.55 + 0.35 * (cx / 1.0)  # ramps from 0.55 to 0.90
        ax.bar(cx, cnt, width=bin_w * 0.88,
               color=shade, alpha=alpha,
               edgecolor="white", linewidth=0.3, zorder=3)

    # Power = 0.80 threshold line
    ax.axvline(0.80, color=CORAL, ls="--", lw=1.8, zorder=5,
               label="Power = 0.80")

    ax.set_xlim(-0.02, 1.02)
    ymax = counts.max() * 1.30
    ax.set_ylim(0, ymax)

    ax.set_xlabel("Per-feature statistical power", fontsize=AXT)
    ax.set_ylabel("Number of tests", fontsize=AXT)

    # Median line
    med = pwr_series.median()
    ax.axvline(med, color="#555555", ls=":", lw=1.2, zorder=4, alpha=0.7)

    # Title
    comp_label = {
        "tissue":   "Tissue comparisons",
        "genotype": "Genotype comparisons",
        "temporal": "Temporal comparisons",
    }.get(label, label.capitalize())
    ax.set_title(comp_label,
                 fontsize=SUB, fontweight="bold", pad=7)

    if subtitle:
        ax.text(0.50, -0.18, subtitle,
                transform=ax.transAxes, fontsize=ANN + 1,
                ha="center", va="top", style="italic", color="#666666")

    ptag(ax, tag, x=tag_x, y=1.06)

# =============================================================================
#  PANEL A — Tissue power (hero: green)
# =============================================================================

power_histogram(
    ax_A, power_series["tissue"],
    label="tissue",
    fill_col=POWER_A_COL,
    tag="A",
    subtitle=None,
)

# =============================================================================
#  PANEL B — Genotype power (grey — explicitly secondary)
# =============================================================================

power_histogram(
    ax_B, power_series["genotype"],
    label="genotype",
    fill_col=POWER_B_COL,
    tag="B",
    subtitle=None,
    tag_x=-0.20,
)

# =============================================================================
#  BOOTSTRAP BOXPLOT HELPER
# =============================================================================

def bootstrap_boxplot(ax, metric_col, ylabel, tag, is_hero=False):
    """
    Draws a bootstrap stability boxplot for one metric.
    metric_col : column name in raw bootstrap file
    ylabel     : y-axis label
    is_hero    : if True, adds G2 modularity Δ annotation
    """
    spine_style(ax, bg="#FDFDFB")

    if metric_col not in boot.columns:
        ax.text(0.5, 0.5, f"Column '{metric_col}'\nnot found",
                transform=ax.transAxes, ha="center", va="center",
                color="red", fontsize=ANN)
        return

    box_data  = []
    colors    = []
    baselines = []

    for i, (tis, gen) in enumerate(GROUPS):
        vals = (boot[(boot["tissue"] == tis) & (boot["genotype"] == gen)]
                [metric_col].dropna().values)
        box_data.append(vals)
        colors.append(STYLE[(tis, gen)])
        bval = get_baseline(tis, gen, metric_col)
        baselines.append(bval)

    # Draw boxplot
    bp = ax.boxplot(
        box_data,
        positions=range(4),
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        zorder=3,
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(color="#888888", linewidth=0.9),
        capprops=dict(color="#888888", linewidth=0.9),
    )
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.60)
        patch.set_linewidth(0.8)

    # Baseline star markers
    for i, bval in enumerate(baselines):
        if not np.isnan(bval):
            ax.plot(i, bval, marker="*", color="black",
                    markersize=12, zorder=6,
                    markeredgecolor="white", markeredgewidth=0.6)

    all_vals = np.concatenate([d for d in box_data if len(d) > 0])
    yspan = all_vals.max() - all_vals.min() if len(all_vals) > 0 else 1
    ymin  = all_vals.min() - yspan * 0.10
    ymax  = all_vals.max() + yspan * 0.20

    ax.set_ylim(ymin, ymax)
    ax.set_xticks(range(4))
    ax.set_xticklabels(GROUP_LABELS, fontsize=TCK)
    ax.set_ylabel(ylabel, fontsize=AXT)

    ptag(ax, tag, x=-0.18, y=1.06)

# =============================================================================
#  PANEL C — Temporal power (grey — secondary, explicitly expected to be low)
# =============================================================================

if ax_C is not None and "temporal" in power_series:
    power_histogram(
        ax_C, power_series["temporal"],
        label="temporal",
        fill_col=POWER_B_COL,
        tag="C",
        subtitle=None,
        tag_x=-0.22,
    )

# =============================================================================
#  PANEL D — Bootstrap density
# =============================================================================

bootstrap_boxplot(ax_D, "density_all",
                  ylabel="Network density", tag="D")
ax_D.set_title("Density", fontsize=SUB, fontweight="bold", pad=6)

# =============================================================================
#  PANEL E — Bootstrap modularity (hero — with Δ annotation)
# =============================================================================

bootstrap_boxplot(ax_E, "modularity_all",
                  ylabel="Louvain modularity Q", tag="E", is_hero=True)
ax_E.set_title("Louvain Modularity", fontsize=SUB, fontweight="bold", pad=6)

# =============================================================================
#  PANEL F — Bootstrap transitivity
# =============================================================================

bootstrap_boxplot(ax_F, "transitivity_lcc",
                  ylabel="Transitivity (LCC)", tag="F")
ax_F.set_title("Transitivity", fontsize=SUB, fontweight="bold", pad=6)

# =============================================================================
#  SHARED LEGENDS
# =============================================================================

# Top row — power threshold legend
leg_top = [
    mpatches.Patch(fc=POWER_A_COL, ec="none", alpha=0.75,
                   label="Tissue comparisons (basis of central claims)"),
    mpatches.Patch(fc=POWER_B_COL, ec="none", alpha=0.75,
                   label="Genotype comparisons (not central claims)"),
    Line2D([0],[0], color=CORAL, ls="--", lw=1.6,
           label="Power = 0.80 threshold"),
    Line2D([0],[0], color="#555555", ls=":", lw=1.2,
           label="Median power"),
]
_leg_ncol = 3 if _has_temporal else 2
_leg_handles_top = leg_top if not _has_temporal else leg_top + [
    mpatches.Patch(fc=POWER_B_COL, ec="none", alpha=0.75,
                   label="Temporal comparisons (not central claims)"),
]
fig.legend(handles=_leg_handles_top,
           loc="upper center", bbox_to_anchor=(0.50, 0.545),
           ncol=_leg_ncol, fontsize=LEG, frameon=True,
           framealpha=0.95, edgecolor=ANNO_EC
           ).get_frame().set_linewidth(0.7)

# Bottom row — bootstrap legend
leg_bot = [
    mpatches.Patch(fc=STYLE[(t, g)], ec="none", alpha=0.65,
                   label=f"{t} {g}")
    for t, g in GROUPS
] + [
    Line2D([0],[0], marker="*", color="black", ms=10,
           lw=0, markeredgecolor="white",
           label="Baseline (all 7 replicates)"),
]
fig.legend(handles=leg_bot,
           loc="lower center", bbox_to_anchor=(0.50, 0.015),
           ncol=5, fontsize=LEG, frameon=True,
           framealpha=0.95, edgecolor=ANNO_EC
           ).get_frame().set_linewidth(0.7)


# =============================================================================
#  SAVE
# =============================================================================

pdf_out = os.path.join(OUTPUT_DIR, "Fig_S9_StatisticalRobustness.pdf")
png_out = os.path.join(OUTPUT_DIR, "Fig_S9_StatisticalRobustness.png")

fig.savefig(pdf_out, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(png_out, dpi=300, bbox_inches="tight", facecolor="white")

print(f"\n{'='*60}")
print(f"  FIGURE S9 COMPLETE")
print(f"  PDF  ->  {pdf_out}")
print(f"  PNG  ->  {png_out}")
print(f"{'='*60}\n")
plt.close(fig)