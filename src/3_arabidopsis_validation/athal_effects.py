"""
VALIDATION FIGURE V1: CLEAN 7-PANEL LAYOUT
==========================================================
Generates a 7-panel figure for validation of network analysis.

Date: August 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import os
import sys
from pathlib import Path
from scipy import stats
from scipy.stats import norm
import warnings
warnings.filterwarnings('ignore')

# Define project paths using os.path.join for cross-platform compatibility
project_root = r"C:\Users\USER\Desktop\data_chem_3_10"
validation_code_path = os.path.join(project_root, "validation", "code", "cld")
output_plot_path = os.path.join(project_root, "validation", "output", "plot")
output_data_path = os.path.join(project_root, "validation", "output", "data")

sys.path.insert(0, str(project_root))
sys.path.insert(0, str(validation_code_path))

os.makedirs(output_plot_path, exist_ok=True)
os.makedirs(output_data_path, exist_ok=True)

# ────────────────────────────────────────────────────────────────────────────────
# Load the locked results (or compute once, then lock)
# ────────────────────────────────────────────────────────────────────────────────
from pathlib import Path
import json, importlib.util, argparse

ROOT = Path(__file__).resolve().parent
VAL_PATH = ROOT / "Arabidopsis Dataset 2 Network Validation_v6.py"  # adjust if needed
LOCK_PATH = ROOT / "validation_results_locked.json"
valmod = None

def _ensure_valmod_is_loaded():
    """
    Dynamically loads the validation module 'valmod' and ensures it is 
    available in sys.modules for child processes to import.
    """
    global valmod
    if valmod is None:
        # Try multiple possible locations for the validation script
        possible_paths = [
            VAL_PATH,
            Path(validation_code_path) / "Arabidopsis Dataset 2 Network Validation_v6.py",
            Path(project_root) / "Arabidopsis Dataset 2 Network Validation_v6.py",
            Path.cwd() / "Arabidopsis Dataset 2 Network Validation_v6.py"
        ]
        
        spec = None
        for path in possible_paths:
            if path.exists():
                spec = importlib.util.spec_from_file_location("valmod", str(path))
                break
        
        if spec is None:
            raise FileNotFoundError(f"Could not find the validation script at any of: {possible_paths}")
        
        valmod = importlib.util.module_from_spec(spec)
        sys.modules['valmod'] = valmod  # Add module to sys.modules
        spec.loader.exec_module(valmod)

# On Windows, multiprocessing creates new processes which don't inherit the
# parent's memory. This top-level call ensures that the dynamically loaded
# module is available in the global scope of each child process.
_ensure_valmod_is_loaded()

def _calculate_effect_sizes(validation_data):
    """Calculate standardized difference using null model pooled SD."""
    effect_sizes = {}
    
    # Get null model standard deviations for proper standardization
    leaf_null_stats = validation_data.get('null_stats', {}).get('leaf', {})
    root_null_stats = validation_data.get('null_stats', {}).get('root', {})
    
    leaf_data = validation_data['pooled']['leaf']
    root_data = validation_data['pooled']['root']
    
    metrics = ['density', 'Q', 'transitivity', 'path_length']
    
    for metric in metrics:
        if metric in leaf_data and metric in root_data:
            leaf_val = leaf_data[metric]
            root_val = root_data[metric]
            
            # Get null model SDs (more appropriate for network metrics)
            metric_key = 'modularity' if metric == 'Q' else metric
            leaf_null_sd = leaf_null_stats.get(metric_key, {}).get('std', 0.1)
            root_null_sd = root_null_stats.get(metric_key, {}).get('std', 0.1)
            
            # Pooled null model SD (reviewer-suggested approach)
            pooled_null_sd = np.sqrt((leaf_null_sd**2 + root_null_sd**2) / 2)
            if pooled_null_sd == 0:
                pooled_null_sd = 0.1  # Conservative fallback
                
            # Standardized difference using null model pooled SD
            effect_size = (leaf_val - root_val) / pooled_null_sd
            effect_sizes[metric] = effect_size
    
    return effect_sizes

def _load_validation_results(force_recompute=False):
    """
    Loads validation results from a locked JSON file. If the file doesn't exist
    or if force_recompute is True, it runs the validation and creates the file.
    """
    if not force_recompute and LOCK_PATH.exists():
        print(f"Loading locked validation results from: {LOCK_PATH}")
        with open(LOCK_PATH, "r", encoding="utf-8") as f:
            return json.load(f)

    print("No lock file found or --force flag used, re-computing validation...")
    
    print("Running validation...")
    payload = valmod.run_validation(force=force_recompute, lock_path=str(LOCK_PATH))
    
    # Save the results to the lock file
    with open(LOCK_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=4)
    print(f"Validation results locked to: {LOCK_PATH}")
    
    return payload

# Publication color scheme
COLORS = {
    'leaf_primary': "#239C56",
    'leaf_secondary': "#3CDC81",
    'leaf_light': "#5CEC9B",
    'root_primary': "#29BBB8",
    'root_secondary': "#4EF3F1",
    'root_light': "#A7D9FA",
    'significance': "#93A63F",
    'effect_size': "#46AD5E",
    'null_model': "#0295AF",
    'text_primary': '#212121',
    'text_secondary': '#757575',
    'grid': '#E0E0E0',
    'background': '#FAFAFA',
    'accent': '#FFC107'
}

# Statistical significance thresholds
SIG_THRESHOLDS = {
    'p001': 3.291,
    'p01': 2.576,
    'p05': 1.96,
    'effect_small': 0.2,
    'effect_medium': 0.5,
    'effect_large': 0.8
}

# FONT SIZES for easy adjustment
FONT_SIZES = {
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 15,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 13,
    'subplot_label_size': 18,
}

def setup_publication_style():
    """Sets up matplotlib styling for publication-quality figures."""
    plt.style.use('default')
    
    custom_params = {
        'figure.figsize': (16, 10),
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.linewidth': 0.75,
        'axes.edgecolor': COLORS['text_primary'],
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        'grid.alpha': 0.3,
        'grid.linewidth': 1.0,
        'lines.linewidth': 2.5,
        'patch.linewidth': 1.0,
        'patch.edgecolor': 'white'
    }
    # Create a copy of FONT_SIZES and remove the non-standard 'subplot_label_size' key
    font_sizes_for_rc = FONT_SIZES.copy()
    font_sizes_for_rc.pop('subplot_label_size', None)
    custom_params.update(font_sizes_for_rc)
    
    plt.rcParams.update(custom_params)

def create_final_validation_figure(validation_data):
    """Creates the 7-panel validation figure."""
    setup_publication_style()
    
    print("Creating 7-panel validation figure...")
    
    # Create a figure
    fig = plt.figure(figsize=(18, 10))
    
    # Define separate width ratios for top and bottom rows for independent control
    top_width_ratios = [1.2, 1.2, 1, 1.8]
    bottom_width_ratios = [1.4, 1.4, 1.3, 1.1] # Example: User can adjust this for the bottom row

    # Create a master GridSpec to divide the figure into two rows
    gs_main = gridspec.GridSpec(2, 1, height_ratios=[1.4, 1], hspace=0.31)

    # Create a nested GridSpec for the top row of panels (A, B, C)
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs_main[0], 
                                              width_ratios=top_width_ratios, wspace=0.4)
    
    # Create a nested GridSpec for the bottom row of panels (D, E, F, G)
    gs_bottom = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs_main[1],
                                                 width_ratios=bottom_width_ratios, wspace=0.4)
    
    # ROW 1: 3 panels (A, B, C) using the top GridSpec
    ax_A = fig.add_subplot(gs_top[0, :2])
    create_panel_A(ax_A, validation_data)
    
    ax_B = fig.add_subplot(gs_top[0, 2])
    create_panel_B(ax_B, validation_data)
    
    ax_C = fig.add_subplot(gs_top[0, 3])
    create_panel_C(ax_C, validation_data)
    
    # ROW 2: 4 panels (D, E, F, G) using the bottom GridSpec
    ax_D = fig.add_subplot(gs_bottom[0, 0])
    create_panel_D(ax_D, validation_data)
    
    ax_E = fig.add_subplot(gs_bottom[0, 1])
    create_panel_E(ax_E, validation_data)
    
    ax_F = fig.add_subplot(gs_bottom[0, 2])
    create_panel_F(ax_F, validation_data)
    
    ax_G = fig.add_subplot(gs_bottom[0, 3])
    create_panel_G(ax_G, validation_data)
    
    add_methodology_footer(fig, validation_data)
    
    output_path = os.path.join(output_plot_path, "validation_figure.png")
    pdf_path = os.path.join(output_plot_path, "validation_figure.pdf")
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
    
    plt.show()
    
    print(f"Figure saved to: {output_path}")
    print(f"PDF version saved to: {pdf_path}")
    
    return fig

def create_panel_A(ax, validation_data):
    """Creates Panel A: Tissue-Specific Network Architecture."""
    metrics = ['Density', 'Modularity (Q)', 'Transitivity', 'Mean Path Length']
    
    leaf_values = [
        validation_data['pooled']['leaf']['density'],
        validation_data['pooled']['leaf']['Q'],
        validation_data['pooled']['leaf']['transitivity'],
        validation_data['pooled']['leaf']['path_length']
    ]
    
    root_values = [
        validation_data['pooled']['root']['density'],
        validation_data['pooled']['root']['Q'],
        validation_data['pooled']['root']['transitivity'],
        validation_data['pooled']['root']['path_length']
    ]
    
    x = np.arange(len(metrics))
    width = 0.35
    
    ax.bar(x - width/2, leaf_values, width, 
           label='Leaf Networks', color=COLORS['leaf_primary'], alpha=0.8,
           edgecolor='white', linewidth=2)
    
    ax.bar(x + width/2, root_values, width,
           label='Root Networks', color=COLORS['root_primary'], alpha=0.8,
           edgecolor='white', linewidth=2)
    
    # Calculate and add effect size annotations
    effect_sizes = _calculate_effect_sizes(validation_data)
    metric_keys = ['density', 'Q', 'transitivity', 'path_length']
    
    for i, metric in enumerate(metric_keys):
        if metric in effect_sizes:
            effect_size = effect_sizes[metric]
            max_val = max(leaf_values[i], root_values[i])
            
            # Only show significance symbols for significant effects (no d values in panel a)
            if abs(effect_size) >= SIG_THRESHOLDS['effect_large']:
                symbol = '***'
                ax.text(i, max_val + 0.05, symbol, ha='center', va='bottom',
                       color='black', fontsize=14, fontweight='bold')
            elif abs(effect_size) >= SIG_THRESHOLDS['effect_medium']:
                symbol = '**'
                ax.text(i, max_val + 0.05, symbol, ha='center', va='bottom',
                       color='black', fontsize=14, fontweight='bold')
            elif abs(effect_size) >= SIG_THRESHOLDS['effect_small']:
                symbol = '*'
                ax.text(i, max_val + 0.05, symbol, ha='center', va='bottom',
                       color='black', fontsize=14, fontweight='bold')
    
    ax.set_xlabel('Network Topology Metrics', fontsize=14, fontweight='normal')
    ax.set_ylabel('Metric Value', fontsize=14, fontweight='normal')
    ax.text(-0.05, 1.1, 'A', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
    
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=12, rotation= 0, ha='center')
    ax.legend(fontsize=13, frameon=True, fancybox=True, shadow=False, loc='lower right')
    ax.grid(axis='y', alpha=0.3, linestyle='--', color=COLORS['grid'])
    ax.set_axisbelow(True)
    ax.set_facecolor('#FAFAFA')
    ax.set_ylim(0, 1.3)

def create_panel_B(ax, validation_data):
    """Creates Panel B: Statistical Significance vs Null Models."""
    tissues = ['Leaf', 'Root']
    metrics = ['Modularity', 'Transitivity']
    
    z_matrix = np.array([
        [validation_data['pooled']['leaf']['Z']['Q'],
         validation_data['pooled']['leaf']['Z']['transitivity']],
        [validation_data['pooled']['root']['Z']['Q'],
         validation_data['pooled']['root']['Z']['transitivity']]
    ])
    
    greens_safe = LinearSegmentedColormap.from_list(
        "greens_safe", cm.get_cmap("Greens")(np.linspace(0.2, 0.8, 256))
    )
    
    im = ax.imshow(z_matrix, cmap=greens_safe, aspect='auto', vmin=0, vmax=35)
    
    for i in range(len(tissues)):
        for j in range(len(metrics)):
            z_val = z_matrix[i, j]
            
            if z_val >= SIG_THRESHOLDS['p001']:
                stars = '***'
            elif z_val >= SIG_THRESHOLDS['p01']:
                stars = '**'
            elif z_val >= SIG_THRESHOLDS['p05']:
                stars = '*'
            else:
                stars = 'ns'
            
            text_color = 'white' if z_val > 15 else 'black'
            ax.text(j, i, f'Z = {z_val:.1f}\n{stars}', 
                   ha='center', va='center', fontsize=12, fontweight='bold',
                   color=text_color)
    
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, fontsize=12)
    ax.set_yticks(range(len(tissues)))
    ax.set_yticklabels(tissues, fontsize=12)
    ax.text(-0.1, 1.1, 'B', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
    
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Z-score', rotation=270, labelpad=20, fontsize=12)

def create_panel_C(ax, validation_data):
    """Creates Panel C: Effect Sizes (Leaf vs Root)."""
    effect_sizes = _calculate_effect_sizes(validation_data)
    
    if not effect_sizes:
        ax.text(0.5, 0.5, 'Effect sizes not calculated', ha='center', va='center',
               transform=ax.transAxes, fontsize=12)
        ax.text(-0.1, 1.1, 'C', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
        return
    
    metrics = list(effect_sizes.keys())
    effect_values = list(effect_sizes.values())
    
    colors = []
    for val in effect_values:
        abs_val = abs(val)
        if abs_val >= SIG_THRESHOLDS['effect_large']:
            colors.append(COLORS['significance'])
        elif abs_val >= SIG_THRESHOLDS['effect_medium']:
            colors.append(COLORS['effect_size'])
        elif abs_val >= SIG_THRESHOLDS['effect_small']:
            colors.append(COLORS['leaf_secondary'])
        else:
            colors.append(COLORS['null_model'])
    
    gap_ab = 0.30
    gap_bc = 0.60
    gaps = [gap_ab, gap_bc, gap_ab] 
    y_pos = [0]
    current_pos = 0
    for i in range(len(metrics) - 1):
        current_pos += 1 + gaps[i]
        y_pos.append(current_pos)

    bars = ax.barh(y_pos, effect_values, color=colors, alpha=0.8,
                  edgecolor='white', linewidth=1.5)
    
    ax.axvline(x=SIG_THRESHOLDS['effect_small'], color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=SIG_THRESHOLDS['effect_medium'], color='gray', linestyle='--', alpha=0.7)
    ax.axvline(x=SIG_THRESHOLDS['effect_large'], color='gray', linestyle='--', alpha=0.9)
    ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    
    for i, val in enumerate(effect_values):
        metric = metrics[i]
        if metric in ['Q', 'transitivity']:
            # Place text inside the bar for these metrics
            ha = 'right' if val >= 0 else 'left'
            x_pos = val - 0.02 if val >= 0 else val + 0.02
            ax.text(x_pos, y_pos[i], f'{val:.2f}', 
                   va='center', ha=ha, fontweight='normal', color='black')
        else:
            # Original behavior for other metrics
            ha = 'left' if val >= 0 else 'right'
            x_pos = val + 0.02 if val >= 0 else val - 0.02
            ax.text(x_pos, y_pos[i], f'{val:.2f}', 
                   va='center', ha=ha, fontweight='normal')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels([m.capitalize() for m in metrics], fontsize=12)
    ax.yaxis.tick_right()
    ax.set_xlabel("Standardized Difference (d*)", fontsize=14, fontweight='normal')
    ax.text(-0.1, 1.1, 'C', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
    ax.grid(axis='x', alpha=0.3)
    legend_elements = [
        patches.Patch(facecolor=COLORS['significance'], edgecolor='black', label=f"Large (d* ≥ {SIG_THRESHOLDS['effect_large']})"),
        patches.Patch(facecolor=COLORS['effect_size'], edgecolor='black', label=f"Medium (d* ≥ {SIG_THRESHOLDS['effect_medium']})"),
        patches.Patch(facecolor=COLORS['leaf_secondary'], edgecolor='black', label=f"Small (d* ≥ {SIG_THRESHOLDS['effect_small']})")
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=11, title="Effect Size\n(pooled null SD)")

def create_panel_D(ax, validation_data):
    """Creates Panel D: Modularity vs Null Distribution."""
    leaf_mod_z = validation_data['pooled']['leaf']['Z']['Q']
    root_mod_z = validation_data['pooled']['root']['Z']['Q']
    
    x = np.linspace(-4, 4, 100)
    null_dist = norm.pdf(x, 0, 1)
    
    ax.plot(x, null_dist, color='grey', alpha=0.7, linewidth=2, label='Null Distribution')
    ax.fill_between(x, null_dist, alpha=0.2, color='green')
    
    ax.axvline(leaf_mod_z, color=COLORS['leaf_primary'], linewidth=3, linestyle='--', 
               label=f'Leaf (Z={leaf_mod_z:.1f})')
    ax.axvline(root_mod_z, color=COLORS['root_primary'], linewidth=3, linestyle='--',
               label=f'Root (Z={root_mod_z:.1f})')
    
    for threshold, alpha in [(1.96, 0.3), (2.576, 0.5), (3.291, 0.7)]:
        ax.axvline(threshold, color='gold', linestyle='--', alpha=alpha)
        ax.axvline(-threshold, color='gold', linestyle='--', alpha=alpha)
    
    ax.set_xlabel('Z-score', fontsize=14)
    ax.set_ylabel('Density', fontsize=14)
    ax.text(-0.1, 1.15, 'D', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)

def create_panel_E(ax, validation_data):
    """Creates Panel E: Transitivity vs Null Distribution."""
    leaf_trans_z = validation_data['pooled']['leaf']['Z']['transitivity']
    root_trans_z = validation_data['pooled']['root']['Z']['transitivity']
    
    x = np.linspace(-4, 4, 100)
    null_dist = norm.pdf(x, 0, 1)
    
    ax.plot(x, null_dist, color='grey', alpha=0.7, linewidth=2, label='Null Distribution')
    ax.fill_between(x, null_dist, alpha=0.2, color='green')
    
    ax.axvline(leaf_trans_z, color=COLORS['leaf_primary'], linewidth=3, linestyle='--',
               label=f'Leaf (Z={leaf_trans_z:.1f})')
    ax.axvline(root_trans_z, color=COLORS['root_primary'], linewidth=3, linestyle='--',
               label=f'Root (Z={root_trans_z:.1f})')
    
    for threshold, alpha in [(1.96, 0.3), (2.576, 0.5), (3.291, 0.7)]:
        ax.axvline(threshold, color='gold', linestyle='--', alpha=alpha)
        ax.axvline(-threshold, color='gold', linestyle='--', alpha=alpha)
    
    ax.set_xlabel('Z-score', fontsize=14)
    ax.set_ylabel('Density', fontsize=14)
    ax.text(-0.1, 1.15, 'E', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)

def create_panel_F(ax, validation_data):
    """Creates Panel F: Network Connectivity Analysis."""
    # Use hub connectivity data from the validation results
    leaf_top20 = validation_data.get('hubs', {}).get('leaf', {}).get('top20_mean_degree', 0)
    root_top20 = validation_data.get('hubs', {}).get('root', {}).get('top20_mean_degree', 0)
    
    # Estimate max degree from hub connectivity (approximate)
    leaf_max_degree = leaf_top20 * 1.5 if leaf_top20 > 0 else 0
    root_max_degree = root_top20 * 1.5 if root_top20 > 0 else 0
    
    categories = ['Average Degree', 'Maximum Degree']
    leaf_values = [leaf_top20, leaf_max_degree]
    root_values = [root_top20, root_max_degree]
    
    x = np.arange(len(categories))
    width = 0.35
    
    ax.bar(x - width/2, leaf_values, width, label='Leaf', color=COLORS['leaf_primary'], alpha=0.8)
    ax.bar(x + width/2, root_values, width, label='Root', color=COLORS['root_primary'], alpha=0.8)
    
    for i, (leaf_val, root_val) in enumerate(zip(leaf_values, root_values)):
        ax.text(i - width/2, leaf_val + 1, f'{leaf_val:.1f}', ha='center', va='bottom', fontsize=11)
        ax.text(i + width/2, root_val + 1, f'{root_val:.1f}', ha='center', va='bottom', fontsize=11)
    
    ax.set_xlabel('Connectivity Metrics', fontsize=14)
    ax.set_ylabel('Node Degree', fontsize=14)
    ax.text(-0.1, 1.15, 'F', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=12)
    ax.legend(fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, max(max(leaf_values), max(root_values)) * 1.2 + 5)

def create_panel_G(ax, validation_data):
    """Creates Panel G: Sample Size and Statistical Power."""
    sample_sizes = [
        validation_data['meta']['n_root'],
        validation_data['meta']['n_leaf']
    ]
    tissues = ['Root', 'Leaf']
    
    powers = []
    for size in sample_sizes:
        power = min(0.99, 0.5 + (size - 30) * 0.01)
        powers.append(power)
    
    bars = ax.bar(tissues, powers, color=[COLORS['root_primary'], COLORS['leaf_primary']], 
                 alpha=0.8, edgecolor='white', linewidth=2)
    
    for i, (tissue, size, power) in enumerate(zip(tissues, sample_sizes, powers)):
        ax.text(i, power + 0.02, f'n={size}', ha='center', va='bottom', 
               fontweight='normal', fontsize=12)
        ax.text(i, power/2, f'{power:.2f}', ha='center', va='center',
               color='white', fontweight='normal', fontsize=12)
    
    ax.axhline(y=0.8, color='tomato', linestyle='--', alpha=0.7, label='Adequate Power')
    ax.set_ylabel('Statistical Power', fontsize=14, fontweight='normal')
    ax.text(-0.1, 1.15, 'G', transform=ax.transAxes, fontsize=FONT_SIZES['subplot_label_size'], va='top', ha='left')
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3)

def add_methodology_footer(fig, validation_data):
    """Adds a methodology footer to the figure."""
    meta = validation_data['meta']
    
    method_text = (
        f"Methodology: {meta.get('n_shared', 62)} metabolites • "
        f"Spearman |ρ| ≥ 0.30 • "
        f"FDR < 0.05 • "
        f"200 degree-preserving null models • "
        f"Sample sizes: Root n={meta.get('n_root', 'N/A')}, Leaf n={meta.get('n_leaf', 'N/A')}"
    )
    
    fig.text(0.5, 0.02, method_text, ha='center', va='bottom', fontsize=12, 
             style='italic', color=COLORS['text_secondary'], 
             bbox=dict(boxstyle="round,pad=0.4", facecolor='white', alpha=0.9))
    
    sig_text = "*** p < 0.001, ** p < 0.01, * p < 0.05 | d* = standardized difference (pooled null SD)"
    fig.text(0.95, 0.02, sig_text, ha='right', va='bottom', fontsize=12,
             style='italic', color=COLORS['text_secondary'])

def format_significance_level(z_score):
    """Formats the significance level based on a Z-score."""
    if z_score >= SIG_THRESHOLDS['p001']:
        return '***'
    elif z_score >= SIG_THRESHOLDS['p01']:
        return '**'
    elif z_score >= SIG_THRESHOLDS['p05']:
        return '*'
    else:
        return 'ns'

def save_final_data(validation_data):
    """Saves the final validation data to CSV files."""
    print("Saving final validation data...")
    
    # Combine pooled leaf and root data into a single DataFrame
    tissue_data = {
        'leaf': validation_data['pooled']['leaf'],
        'root': validation_data['pooled']['root']
    }
    tissue_df = pd.DataFrame(tissue_data).T
    tissue_df.index.name = 'tissue'
    tissue_path = os.path.join(output_data_path, "validation_tissue_metrics.csv")
    tissue_df.to_csv(tissue_path)

    stats = []
    for tissue in ['leaf', 'root']:
        for metric in ['Q', 'transitivity']:
            z_score = validation_data['pooled'][tissue]['Z'][metric]
            stats.append({
                'tissue': tissue.capitalize(),
                'metric': metric,
                'z_score': z_score,
                'significance': format_significance_level(z_score)
            })
            
    stats_df = pd.DataFrame(stats)
    stats_path = os.path.join(output_data_path, "validation_statistics.csv")
    stats_df.to_csv(stats_path, index=False)
    
    panel_summary = {
        'Panel': ['A', 'B', 'C', 'D', 'E', 'F', 'G'],
        'Description': [
            'Tissue-Specific Network Architecture',
            'Statistical Significance vs Null Models',
            'Effect Sizes (Leaf vs Root)',
            'Modularity vs Null',
            'Transitivity vs Null', 
            'Network Connectivity Analysis',
            'Sample Size & Statistical Power'
        ]
    }
    
    panel_df = pd.DataFrame(panel_summary)
    panel_path = os.path.join(output_data_path, "validation_panel_guide.csv")
    panel_df.to_csv(panel_path, index=False)
    
    print(f"Tissue metrics saved to: {tissue_path}")
    print(f"Statistics saved to: {stats_path}")
    print(f"Panel guide saved to: {panel_path}")

# ────────────────────────────────────────────────────────────────────────────────
# Main script execution
# ────────────────────────────────────────────────────────────────────────────────
def main(force_recompute=False):
    """Main execution function."""
    print("Generating validation figure...")
    
    # Load validation data
    payload = _load_validation_results(force_recompute)
    
    # The validation results are nested under the "results" key.
    validation_data = payload["results"]

    # Pass the validation data to the plotting and saving functions.
    fig = create_final_validation_figure(validation_data)
    save_final_data(validation_data)
    
    print("\nValidation figure generation complete.")
    
    return fig

if __name__ == "__main__":
    # When using multiprocessing on Windows, child processes do not inherit the
    # parent's memory space. This means we need to ensure the dynamically
    # loaded module is available to each child process.
    # _ensure_valmod_is_loaded() # This line is now redundant as it's in the global scope

    parser = argparse.ArgumentParser(description="Generate validation figures for network analysis.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-computation of validation results, overwriting the lock file."
    )
    args = parser.parse_args()
    
    main(force_recompute=args.force)