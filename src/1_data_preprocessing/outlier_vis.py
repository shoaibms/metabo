"""
Visualizes the impact of outlier imputation on metabolomics data.

This script compares metabolomics data before and after outlier imputation
to assess the magnitude of changes for different metabolite clusters. It
calculates the relative impact of imputation and generates an error bar plot
to visualize the differences for the most affected variables.
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Any

# --- Plotting Constants ---
PRE_IMPUTATION_COLOR = '#65c792'
POST_IMPUTATION_COLOR = '#0e3823'
PRE_IMPUTATION_ECOLOR = '#7fcdb0'
POST_IMPUTATION_ECOLOR = '#2f5240'
FIGURE_SIZE = (12, 8)
TITLE_FONT_SIZE = 18
AXIS_LABEL_FONT_SIZE = 16
TICK_LABEL_FONT_SIZE = 12
LEGEND_FONT_SIZE = 14

# --- Configuration ---
CONFIG = {
    'pre_imputation_path': Path(r'C:\Users\ms\Desktop\data_chem\data\impute_missing_value\p_l_rf.csv'),
    'post_imputation_path': Path(r'C:\Users\ms\Desktop\data_chem\data\outlier_fixed\p_l_if.csv'),
    'cluster_prefix': 'P_Cluster',
    'n_top': 30
}

def load_datasets(pre_path: Path, post_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Loads pre- and post-imputation datasets from CSV files.

    Args:
        pre_path: The file path for the pre-imputation data.
        post_path: The file path for the post-imputation data.

    Returns:
        A tuple containing the pre- and post-imputation DataFrames.
    """
    data_before = pd.read_csv(pre_path)
    data_after = pd.read_csv(post_path)
    return data_before, data_after


def get_cluster_columns(data: pd.DataFrame, prefix: str) -> List[str]:
    """Extracts column names that match a given prefix.

    Args:
        data: The DataFrame to search for columns.
        prefix: The prefix to identify cluster columns.

    Returns:
        A list of column names that start with the specified prefix.
    """
    return [col for col in data.columns if col.startswith(prefix)]


def calculate_impact_metrics(
    data_before: pd.DataFrame, data_after: pd.DataFrame, cluster_cols: List[str]
) -> pd.DataFrame:
    """Calculates imputation impact metrics for metabolite clusters.

    Args:
        data_before: DataFrame with data before outlier imputation.
        data_after: DataFrame with data after outlier imputation.
        cluster_cols: A list of column names for the metabolite clusters to analyze.

    Returns:
        A DataFrame with statistics and imputation impact, sorted by absolute impact.
    """
    before_stats = data_before[cluster_cols].agg(["mean", "std"]).T
    before_stats.columns = ["Before", "Before_Std"]

    after_stats = data_after[cluster_cols].agg(["mean", "std"]).T
    after_stats.columns = ["After", "After_Std"]

    impact_df = pd.concat([before_stats, after_stats], axis=1)

    # Calculate relative impact, handling division by zero.
    impact = (impact_df["After"] - impact_df["Before"]) / impact_df["Before"]
    impact_df["Impact"] = impact.replace([np.inf, -np.inf], 0).fillna(0)
    impact_df["AbsImpact"] = impact_df["Impact"].abs()

    impact_df = impact_df.sort_values("AbsImpact", ascending=False)
    impact_df = impact_df.reset_index().rename(columns={"index": "Variable"})
    
    return impact_df


def plot_impact_comparison(impact_df: pd.DataFrame, n_top: int):
    """Creates an error bar plot to compare pre- and post-imputation values.

    Args:
        impact_df: DataFrame containing the impact metrics.
        n_top: The number of top impacted variables to plot.
    """
    top_impact = impact_df.head(n_top)

    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=FIGURE_SIZE)

    # Plot pre- and post-imputation values with error bars
    plt.errorbar(
        range(n_top), top_impact['Before'], yerr=top_impact['Before_Std'],
        fmt='o', color=PRE_IMPUTATION_COLOR, ecolor=PRE_IMPUTATION_ECOLOR,
        capsize=5, label='Pre-Imputation', alpha=0.8
    )
    plt.errorbar(
        range(n_top), top_impact['After'], yerr=top_impact['After_Std'],
        fmt='o', color=POST_IMPUTATION_COLOR, ecolor=POST_IMPUTATION_ECOLOR,
        capsize=5, label='Post-Imputation', alpha=0.8
    )

    # Add reference line and styling
    plt.axhline(y=0, color='gray', linestyle='--')
    plt.ylabel('Mean Value', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.xlabel('Metabolite Clusters', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.title(f'Top {n_top} Most Impacted Variables from Outlier Imputation', fontsize=TITLE_FONT_SIZE)
    plt.xticks(range(n_top), top_impact['Variable'], rotation=90, fontsize=TICK_LABEL_FONT_SIZE)
    plt.yticks(fontsize=TICK_LABEL_FONT_SIZE)
    plt.legend(fontsize=LEGEND_FONT_SIZE)

    plt.tight_layout()
    plt.show()


def main(config: Dict[str, Any]):
    """Main function to run the outlier imputation impact analysis."""
    # Load and process data
    data_before, data_after = load_datasets(config['pre_imputation_path'], config['post_imputation_path'])
    cluster_cols = get_cluster_columns(data_before, config['cluster_prefix'])
    impact_df = calculate_impact_metrics(data_before, data_after, cluster_cols)

    # --- Visualization ---
    # Generate and display the plot
    plot_impact_comparison(impact_df, config['n_top'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize the impact of outlier imputation.")
    parser.add_argument("--pre_path", type=Path, default=CONFIG['pre_imputation_path'],
                        help="Path to the pre-imputation data file.")
    parser.add_argument("--post_path", type=Path, default=CONFIG['post_imputation_path'],
                        help="Path to the post-imputation data file.")
    parser.add_argument("--prefix", type=str, default=CONFIG['cluster_prefix'],
                        help="Prefix for cluster columns.")
    parser.add_argument("--n_top", type=int, default=CONFIG['n_top'],
                        help="Number of top impacted variables to plot.")
    args = parser.parse_args()

    current_config = {
        'pre_imputation_path': args.pre_path,
        'post_imputation_path': args.post_path,
        'cluster_prefix': args.prefix,
        'n_top': args.n_top
    }
    main(current_config)
