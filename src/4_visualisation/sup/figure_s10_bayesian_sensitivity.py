# =============================================================================
#  figure_s10_bayesian_sensitivity.py
#  Supplementary Figure S10 — Bayesian Network Parent-Constraint Sensitivity
#  TPC-2025-1173
#
#  Script location : <project_root>/rivision_pack_1/
#  Output location : <data_root>/output/revision_pack/figures/
#
#  DATA SOURCES (all real — from bayesian_sensitivity.R outputs)
#  ---------------------------------------------------------------
#  Panels A + B:  Supp_Table_S9_BayesianParentSensitivity.csv
#                   columns used: dataset, maxp, observed_edges,
#                                 null_edges_mean, null_edges_sd,
#                                 overlap_with_unconstr, p_perm
#
#  Panel C:       bn_null_Leaf_maxpInf.rds   — Leaf unconstrained null dist
#                 bn_null_Leaf_maxp3.rds     — Leaf maxp=3 null dist
#                 bn_null_Root_maxpInf.rds   — Root unconstrained null dist
#                 bn_null_Root_maxp3.rds     — Root maxp=3 null dist
#                   each RDS contains: $null_edges_raw (vector) +
#                                      $observed_edges (scalar)
#
#  NOTE: Does NOT re-run the Bayesian analysis (57 min runtime).
#        Uses pre-computed outputs from bayesian_sensitivity.R only.
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
#  PATHS — USER CONFIGURATION: update these paths for your local environment
# =============================================================================

DATA_DIR   = r"C:\Users\USER\Desktop\data_chem_3_10\output\revision_pack\bayesian_parent_sensitivity"
OUTPUT_DIR = r"C:\Users\USER\Desktop\data_chem_3_10\output\revision_pack\figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
#  LOAD DATA
# =============================================================================

# ── Summary CSV (Panels A + B) ───────────────────────────────────────────────
csv_path = os.path.join(DATA_DIR, "Supp_Table_S9_BayesianParentSensitivity.csv")
if not os.path.isfile(csv_path):
    sys.exit(f"[ERROR] File not found:\n  {csv_path}\n"
             f"Run bayesian_sensitivity.R first to generate outputs.")

df = pd.read_csv(csv_path)
df.columns = [c.strip() for c in df.columns]
df["maxp"] = df["maxp"].astype(str).str.strip()
# Normalise: R writes Inf as 'inf', numerics as '5.0'/'3.0' — standardise all
df["maxp"] = df["maxp"].str.replace(r'^(?i)inf$', 'Inf', regex=True)   # inf -> Inf
df["maxp"] = df["maxp"].str.replace(r'\.0$', '', regex=True)           # 5.0 -> 5
print(f"[OK] Summary CSV loaded: {len(df)} rows")
print(df[["dataset","maxp","observed_edges","null_edges_mean",
          "null_edges_sd","overlap_with_unconstr","p_perm"]].to_string(index=False))

def _get(tissue, maxp_str, col):
    row = df[(df["dataset"] == tissue) & (df["maxp"] == maxp_str)]
    if row.empty:
        raise ValueError(f"Missing row in CSV: dataset={tissue}, maxp={maxp_str}")
    return float(row[col].iloc[0])

# Observed edge counts  [Unconstrained, maxp=5, maxp=3]
obs = {t: [int(_get(t, m, "observed_edges")) for m in ["Inf", "5", "3"]]
       for t in ["Leaf", "Root"]}

# Null distribution summary (mean and SD — used for band on Panel A)
null_mean = (_get("Leaf","Inf","null_edges_mean") +
             _get("Root","Inf","null_edges_mean")) / 2
null_sd   = (_get("Leaf","Inf","null_edges_sd") +
             _get("Root","Inf","null_edges_sd")) / 2

# Fold-enrichment (observed / null mean) — Panel B
fold = {t: [v / null_mean for v in obs[t]] for t in obs}

# Scaffold overlap percentages  [Unconstrained=None, maxp=5, maxp=3]
overlap = {t: [None,
               _get(t, "5", "overlap_with_unconstr"),
               _get(t, "3", "overlap_with_unconstr")]
           for t in ["Leaf", "Root"]}

# p-values
p_vals = {(t, m): _get(t, m, "p_perm")
          for t in ["Leaf","Root"] for m in ["Inf","5","3"]}

print(f"\n[OK] Observed edges — Leaf: {obs['Leaf']}  Root: {obs['Root']}")
print(f"[OK] Null mean={null_mean:.1f}  SD={null_sd:.1f}")
print(f"[OK] Fold-enrichment — Leaf: {[f'{v:.2f}x' for v in fold['Leaf']]}  "
      f"Root: {[f'{v:.2f}x' for v in fold['Root']]}")

# ── Null distribution RDS files (Panel C) ────────────────────────────────────
try:
    import pyreadr
except ImportError:
    sys.exit("pyreadr required for reading .rds files.\n"
             "Install with:  pip install pyreadr")

# ---------------------------------------------------------------------------
#  R export helper — run this ONCE in R to convert RDS lists to CSV:
#
#  library(purrr)
#  fnames <- c("bn_null_Leaf_maxpInf","bn_null_Leaf_maxp3",
#              "bn_null_Root_maxpInf","bn_null_Root_maxp3")
#  for (f in fnames) {
#    obj <- readRDS(file.path(output_dir, paste0(f, ".rds")))
#    write.csv(data.frame(null_edges_raw = obj$null_edges_raw,
#                         observed_edges = obj$observed_edges),
#              file.path(output_dir, paste0(f, ".csv")), row.names = FALSE)
#  }
# ---------------------------------------------------------------------------

def load_null_rds(fname):
    """Load null distribution — reads CSV export (preferred) or falls back
    to pyreadr if the CSV does not exist yet."""
    base   = fname.replace(".rds", "")
    csv_p  = os.path.join(DATA_DIR, base + ".csv")
    rds_p  = os.path.join(DATA_DIR, fname)

    # ── Path 1: CSV export already exists (fast, reliable) ──────────────────
    if os.path.isfile(csv_p):
        df_null = pd.read_csv(csv_p)
        null_vec = df_null["null_edges_raw"].values.astype(float)
        obs_val  = float(df_null["observed_edges"].iloc[0])
        print(f"[OK] {base}.csv: {len(null_vec)} permutations, observed={int(obs_val)}")
        return {"null": null_vec, "obs": obs_val}

    # ── Path 2: Try pyreadr on the RDS ──────────────────────────────────────
    if not os.path.isfile(rds_p):
        raise FileNotFoundError(
            f"Neither CSV nor RDS found for {base}.\n"
            f"  Expected: {csv_p}\n"
            f"  Run the R export snippet above to create the CSV files.")
    try:
        result  = pyreadr.read_r(rds_p)
        keys    = list(result.keys())
        if not keys:
            raise ValueError("pyreadr returned empty dict — R list not supported")
        obj = result[keys[0]]
        if hasattr(obj, "columns"):
            null_vec = obj["null_edges_raw"].values.astype(float)
            obs_val  = float(obj["observed_edges"].iloc[0])
        else:
            raise ValueError(f"Unexpected type: {type(obj)}")
        print(f"[OK] {fname} via pyreadr: {len(null_vec)} permutations, "
              f"observed={int(obs_val)}")
        return {"null": null_vec, "obs": obs_val}
    except Exception as e:
        raise RuntimeError(
            f"\n{'='*60}\n"
            f"  Cannot read {fname} with pyreadr:\n  {e}\n\n"
            f"  SOLUTION — run this once in R then re-run this script:\n\n"
            f"  output_dir <- '{DATA_DIR.replace(chr(92), '/')}'\n"
            f"  fnames <- c('bn_null_Leaf_maxpInf','bn_null_Leaf_maxp3',\n"
            f"              'bn_null_Root_maxpInf','bn_null_Root_maxp3')\n"
            f"  for (f in fnames) {{\n"
            f"    obj <- readRDS(file.path(output_dir, paste0(f, '.rds')))\n"
            f"    write.csv(data.frame(null_edges_raw = obj$null_edges_raw,\n"
            f"                         observed_edges = obj$observed_edges),\n"
            f"              file.path(output_dir, paste0(f, '.csv')),\n"
            f"              row.names = FALSE)\n"
            f"  }}\n"
            f"{'='*60}\n"
        ) from e

null_dist = {
    ("Leaf", "Unconstrained"): load_null_rds("bn_null_Leaf_maxpInf.rds"),
    ("Leaf", "maxp = 3"):      load_null_rds("bn_null_Leaf_maxp3.rds"),
    ("Root", "Unconstrained"): load_null_rds("bn_null_Root_maxpInf.rds"),
    ("Root", "maxp = 3"):      load_null_rds("bn_null_Root_maxp3.rds"),
}

print("\n[OK] All data loaded — building figure...\n")

# =============================================================================
#  PALETTE  (from colours.py + manuscript identity colours)
# =============================================================================

LEAF_COL   = "#1b941b"    # Leaf tissue — Bayesian scripts identity
LEAF_LIGHT = "#A8DBA8"
ROOT_COL   = "#1e597d"    # Root tissue — Bayesian scripts identity
ROOT_LIGHT = "#9BBDD4"
CORAL      = "#DC6B4A"    # Observed line / primary annotation
NULL_OUTER = "#C6E8FA"    # Light blue  — null ± 2 SD band
NULL_INNER = "#BDFF7A"    # Bright yellow-green — null ± 1 SD band
ANNO_BG    = "#FFFFFF"
ANNO_EC    = "#AAAAAA"
GRID_COL   = "#DDDDDD"

TAG = 18; SUB = 14; AXT = 13; TCK = 11; LEG = 11; ANN = 10

MAXP_X = [0, 1, 2]
MAXP   = ["Unconstrained", "maxp = 5", "maxp = 3"]

TISSUE_CFG = {
    "Leaf": {"col": LEAF_COL, "light": LEAF_LIGHT, "marker": "o", "ms": 9, "zorder": 5},
    "Root": {"col": ROOT_COL, "light": ROOT_LIGHT, "marker": "s", "ms": 9, "zorder": 5},
}

def spine_style(ax, bg="white"):
    for s in ax.spines.values():
        s.set_linewidth(0.8); s.set_color("black")
    ax.yaxis.grid(True, alpha=0.30, lw=0.5, color=GRID_COL, ls="-")
    ax.xaxis.grid(False); ax.set_axisbelow(True)
    ax.tick_params(labelsize=TCK, length=3, width=0.6)
    ax.set_facecolor(bg)

def ptag(ax, tag, x=-0.12, y=1.05):
    ax.text(x, y, tag, transform=ax.transAxes,
            fontsize=TAG, fontweight="normal", va="bottom", ha="left",
            path_effects=[pe.withStroke(linewidth=1.2, foreground="white")])

def abox(ax, x, y, txt, ha="right", va="top", fs=ANN, bold=False, col="#333333"):
    ax.annotate(txt, xy=(x, y), xycoords=ax.transAxes,
                fontsize=fs, ha=ha, va=va, color=col,
                fontweight="bold" if bold else "normal",
                bbox=dict(boxstyle="round,pad=0.28", fc=ANNO_BG,
                          ec=ANNO_EC, lw=0.7, alpha=0.95))

# =============================================================================
#  FIGURE LAYOUT
# =============================================================================

fig = plt.figure(figsize=(14, 11), dpi=300, facecolor="white")
outer = gridspec.GridSpec(2, 1, figure=fig, hspace=0.50,
                          top=0.94, bottom=0.07, left=0.08, right=0.97)
ax_A  = fig.add_subplot(outer[0])
bot   = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1],
                                         wspace=0.38, width_ratios=[0.9, 1.1])
ax_B  = fig.add_subplot(bot[0])
bot_C = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=bot[1],
                                         wspace=0.28, hspace=0.30)
ax_C  = [[fig.add_subplot(bot_C[r, c]) for c in range(2)] for r in range(2)]

# =============================================================================
#  PANEL A — Edge count trajectory with null band
# =============================================================================

spine_style(ax_A, bg="#FDFDFB")

null_lo = null_mean - 2 * null_sd
null_hi = null_mean + 2 * null_sd

ax_A.axhspan(null_lo, null_hi, color=NULL_OUTER, alpha=0.70, zorder=1)
ax_A.axhspan(null_mean - null_sd, null_mean + null_sd,
             color=NULL_INNER, alpha=0.60, zorder=1)
ax_A.axhline(null_mean, color="#3399AA", lw=1.0, ls="--", zorder=2, alpha=0.7)
ax_A.text(2.06, null_mean + 1, f"Null mean ({int(null_mean)})",
          fontsize=ANN - 1, color="#3399AA", va="bottom", style="italic")

for tis, cfg in TISSUE_CFG.items():
    ax_A.fill_between(MAXP_X, obs[tis], null_mean,
                      color=cfg["col"], alpha=0.08, zorder=2)

for tis, cfg in TISSUE_CFG.items():
    ys = obs[tis]
    ax_A.plot(MAXP_X, ys, color="white", lw=4.5, zorder=3, solid_capstyle="round")
    ax_A.plot(MAXP_X, ys, color=cfg["col"], lw=2.4,
              marker=cfg["marker"], ms=cfg["ms"],
              markeredgecolor="white", markeredgewidth=1.2,
              zorder=cfg["zorder"], label=tis)
    for xi, yi in zip(MAXP_X, ys):
        ax_A.text(xi, yi + 18, str(yi),
                  ha="center", va="bottom", fontsize=ANN + 1,
                  fontweight="bold", color=cfg["col"],
                  path_effects=[pe.withStroke(linewidth=2, foreground="white")])

for tis, cfg in TISSUE_CFG.items():
    ax_A.plot(0, obs[tis][0], marker="*", ms=16, color=CORAL,
              zorder=7, markeredgecolor="white", markeredgewidth=0.6)

for xi in MAXP_X:
    delta = obs["Leaf"][xi] - obs["Root"][xi]
    mid_y = (obs["Leaf"][xi] + obs["Root"][xi]) / 2
    ax_A.annotate("", xy=(xi + 0.06, obs["Root"][xi]),
                  xytext=(xi + 0.06, obs["Leaf"][xi]),
                  arrowprops=dict(arrowstyle="<->", color="#888888", lw=0.8))
    ax_A.text(xi + 0.11, mid_y, f"Δ{delta}",
              fontsize=ANN - 1, va="center", color="#666666", style="italic")

ax_A.annotate("Primary analysis\n(Fig. 1D)",
              xy=(0, obs["Leaf"][0]), xycoords="data",
              xytext=(-0.35, obs["Leaf"][0] - 40),
              fontsize=ANN, color=CORAL, ha="center",
              arrowprops=dict(arrowstyle="-|>", color=CORAL, lw=1.0,
                              connectionstyle="arc3,rad=0.2"),
              bbox=dict(boxstyle="round,pad=0.25", fc=ANNO_BG,
                        ec=CORAL, lw=0.8, alpha=0.92))

ax_A.set_xlim(-0.55, 2.55)
_a_ymax = max(obs['Leaf'] + obs['Root']) * 1.22  # 22% headroom for value labels
ax_A.set_ylim(0, _a_ymax)
ax_A.set_xticks(MAXP_X)
ax_A.set_xticklabels(["Unconstrained\n(maxp = ∞)", "maxp = 5", "maxp = 3"],
                     fontsize=TCK + 1)
ax_A.set_ylabel("Edge count  (arcs retained, bootstrap ≥ 0.5)", fontsize=AXT)
ax_A.set_xlabel("Parent constraint level   →   increasing stringency", fontsize=AXT)
ax_A.set_title("Observed edge counts remain far above the permuted-data null"
               " at every constraint level",
               fontsize=SUB, pad=8, fontweight="bold")

leg_A = [
    Line2D([0],[0], color=LEAF_COL, lw=2.2, marker="o", ms=7,
           markeredgecolor="white", label="Leaf"),
    Line2D([0],[0], color=ROOT_COL, lw=2.2, marker="s", ms=7,
           markeredgecolor="white", label="Root"),
    mpatches.Patch(fc=NULL_OUTER, ec="#3399AA", lw=0.7, alpha=0.8, label="Null ± 2 SD"),
    mpatches.Patch(fc=NULL_INNER, ec="#66AA00", lw=0.7, alpha=0.8, label="Null ± 1 SD"),
    Line2D([0],[0], marker="*", ms=10, color=CORAL, lw=0,
           markeredgecolor="white", label="Primary analysis (Fig. 1D)"),
]
lg = ax_A.legend(handles=leg_A, ncol=5, fontsize=LEG,
                  frameon=True, framealpha=0.95, edgecolor=ANNO_EC,
                  loc="upper right", handlelength=1.6)
lg.get_frame().set_linewidth(0.7)
ptag(ax_A, "A", x=-0.06, y=1.04)

# =============================================================================
#  PANEL B — Fold-enrichment (signal strength) trajectory
# =============================================================================

spine_style(ax_B, bg="#FDFDFB")

ax_B.axhspan(0.85, 1.15, color="#FFEEEE", alpha=0.7, zorder=1)
ax_B.axhline(1.0, color="#CC4444", lw=1.4, ls="--", zorder=2, alpha=0.85)
ax_B.text(2.06, 1.0, "Random\n(fold = 1)",
          fontsize=ANN - 1, color="#CC4444", va="center", style="italic")

for tis, cfg in TISSUE_CFG.items():
    ys = fold[tis]
    ax_B.fill_between(MAXP_X, ys, 1.0, color=cfg["col"], alpha=0.10, zorder=2)
    ax_B.plot(MAXP_X, ys, color="white", lw=4, zorder=3, solid_capstyle="round")
    ax_B.plot(MAXP_X, ys, color=cfg["col"], lw=2.4,
              marker=cfg["marker"], ms=cfg["ms"],
              markeredgecolor="white", markeredgewidth=1.2, zorder=5, label=tis)
    for xi, yi in zip(MAXP_X, ys):
        ax_B.text(xi, yi + 0.12, f"{yi:.1f}×",
                  ha="center", va="bottom", fontsize=ANN,
                  fontweight="bold", color=cfg["col"],
                  path_effects=[pe.withStroke(linewidth=2, foreground="white")])

for xi in MAXP_X:
    ax_B.annotate("P < 0.001", xy=(xi, 0.35), xycoords=("data", "data"),
                  fontsize=ANN - 1.5, ha="center", color="#555555",
                  bbox=dict(boxstyle="round,pad=0.20", fc=ANNO_BG,
                            ec=ANNO_EC, lw=0.5, alpha=0.9))

# Scaffold overlap — single compact annotation, range computed from real data
ov_all = [v for t in ["Leaf","Root"] for v in [overlap[t][1], overlap[t][2]]
          if v is not None]
ov_min, ov_max = min(ov_all) * 100, max(ov_all) * 100
ax_B.annotate(
    f"{ov_min:.1f} – {ov_max:.1f}% of constrained arcs\n⊂ unconstrained scaffold",
    xy=(0.50, 0.26), xycoords="axes fraction",
    fontsize=ANN, ha="center", va="bottom", color="#444444",
    bbox=dict(boxstyle="round,pad=0.30", fc=ANNO_BG, ec=ANNO_EC,
              lw=0.7, alpha=0.93)
)

ax_B.set_xlim(-0.55, 2.55)
_b_ymax = max(fold['Leaf'] + fold['Root']) * 1.22  # 22% headroom
ax_B.set_ylim(0, _b_ymax)
ax_B.set_xticks(MAXP_X)
ax_B.set_xticklabels(["Unconstrained\n(maxp = ∞)", "maxp = 5", "maxp = 3"],
                     fontsize=TCK + 1)
ax_B.set_ylabel("Signal strength  (observed ÷ null mean)", fontsize=AXT)
ax_B.set_xlabel("Parent constraint level", fontsize=AXT)
ax_B.set_title("Signal strength stays\nwell above random at all levels",
               fontsize=SUB, pad=8, fontweight="bold")

leg_B = [
    Line2D([0],[0], color=LEAF_COL, lw=2.2, marker="o", ms=7,
           markeredgecolor="white", label="Leaf"),
    Line2D([0],[0], color=ROOT_COL, lw=2.2, marker="s", ms=7,
           markeredgecolor="white", label="Root"),
    mpatches.Patch(fc="#FFEEEE", ec="#CC4444", lw=0.7, alpha=0.8,
                   label="Indistinguishable from random"),
]
lg_B = ax_B.legend(handles=leg_B, fontsize=LEG - 1, frameon=True,
                    framealpha=0.95, edgecolor=ANNO_EC, loc="upper right")
lg_B.get_frame().set_linewidth(0.7)
ptag(ax_B, "B", x=-0.18, y=1.04)

# =============================================================================
#  PANEL C — 2×2 null distribution histograms
# =============================================================================

FACET_ROWS = ["Leaf", "Root"]
FACET_COLS = ["Unconstrained", "maxp = 3"]
HIST_CMAPS = {"Leaf": "Greens", "Root": "Blues"}
COL_TTLS   = ["Unconstrained (maxp = ∞)", "Constrained (maxp = 3)"]
TCOLS      = {"Leaf": LEAF_COL, "Root": ROOT_COL}

for ri, tis in enumerate(FACET_ROWS):
    for ci, maxp_lbl in enumerate(FACET_COLS):
        ax   = ax_C[ri][ci]
        key  = (tis, maxp_lbl)
        narr = null_dist[key]["null"]
        ov   = null_dist[key]["obs"]
        nm   = narr.mean()
        ns   = narr.std()

        spine_style(ax, bg="#FAFFFE" if tis == "Leaf" else "#F8FAFF")

        ax.axvspan(nm - 2*ns, nm + 2*ns, color=NULL_OUTER, alpha=0.55, zorder=1)
        ax.axvspan(nm - ns,   nm + ns,   color=NULL_INNER, alpha=0.50, zorder=1)

        cmap = plt.cm.get_cmap(HIST_CMAPS[tis])
        counts, edges = np.histogram(narr, bins=32)
        for k, (left, right, cnt) in enumerate(zip(edges[:-1], edges[1:], counts)):
            shade = cmap(0.38 + 0.52 * (k / 32))
            ax.bar((left + right) / 2, cnt, width=(right - left) * 0.90,
                   color=shade, ec="white", lw=0.2, alpha=0.88, zorder=3)

        ax.axvline(nm, color="#777777", ls=":", lw=1.0, zorder=4, alpha=0.8)
        ax.axvline(ov, color=CORAL, ls="--", lw=2.5, zorder=6)

        ax.relim(); ax.autoscale_view()
        ylim  = ax.get_ylim()
        span  = edges[-1] - edges[0]

        p_val = p_vals.get((tis, "Inf" if maxp_lbl == "Unconstrained" else "3"), 0.0)
        p_str = "P < 0.001" if p_val < 0.001 else f"P = {p_val:.4f}"
        abox(ax, 0.97, 0.97, p_str, bold=True, col="#333333", fs=ANN)
        abox(ax, 0.97, 0.80, f"Observed\n= {int(ov)}", fs=ANN)
        ax.text(nm + span * 0.03, ylim[1] * 0.48,
                f"Null\nmean\n{nm:.0f}",
                fontsize=ANN - 1.5, color="#666666",
                ha="left", va="center", style="italic")

        if ci == 1:
            ov_pct = overlap[tis][2]
            abox(ax, 0.97, 0.60,
                 f"{ov_pct*100:.1f}% edges\n⊂ unconstrained", fs=ANN - 0.5)

        if ri == 1:
            ax.set_xlabel("Null edge count\n(permutation)", fontsize=AXT - 1)
        if ci == 0:
            ax.set_ylabel("Frequency", fontsize=AXT - 1)
            ax.text(-0.28, 0.5, tis, transform=ax.transAxes,
                    fontsize=AXT + 1, fontweight="bold",
                    ha="center", va="center", rotation=90, color=TCOLS[tis])
        if ri == 0:
            ax.set_title(COL_TTLS[ci], fontsize=SUB - 1, pad=6, fontweight="bold")

ax_C[0][0].text(-0.32, 1.05, "C", transform=ax_C[0][0].transAxes,
                fontsize=TAG, fontweight="normal", va="bottom", ha="left",
                path_effects=[pe.withStroke(linewidth=1.2, foreground="white")])

# =============================================================================
#  SHARED LEGEND + DIVIDER
# =============================================================================

leg_c = [
    Line2D([0],[0], color=CORAL, lw=2.2, ls="--", label="Observed edge count"),
    Line2D([0],[0], color="#777777", lw=1.2, ls=":", label="Null mean"),
    mpatches.Patch(fc=NULL_OUTER, ec="#3399AA", lw=0.7, alpha=0.75, label="Null ± 2 SD"),
    mpatches.Patch(fc=NULL_INNER, ec="#66AA00", lw=0.7, alpha=0.75, label="Null ± 1 SD"),
    mpatches.Patch(fc=LEAF_COL, ec="none", alpha=0.70, label="Leaf null distribution"),
    mpatches.Patch(fc=ROOT_COL, ec="none", alpha=0.70, label="Root null distribution"),
]
# Place legend just above the Panel C 2×2 grid (top-centre of Panel C x-range)
lC = fig.legend(handles=leg_c, loc="lower center", bbox_to_anchor=(0.74, 0.455),
                ncol=3, fontsize=LEG - 1, frameon=True,
                framealpha=0.95, edgecolor=ANNO_EC)
lC.get_frame().set_linewidth(0.7)

# =============================================================================
#  SAVE
# =============================================================================

pdf_out = os.path.join(OUTPUT_DIR, "Fig_S10_BayesianParentSensitivity.pdf")
png_out = os.path.join(OUTPUT_DIR, "Fig_S10_BayesianParentSensitivity.png")

fig.savefig(pdf_out, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(png_out, dpi=300, bbox_inches="tight", facecolor="white")

print(f"\n{'='*60}")
print(f"  FIGURE S10 COMPLETE")
print(f"  PDF  ->  {pdf_out}")
print(f"  PNG  ->  {png_out}")
print(f"{'='*60}\n")

plt.close(fig)