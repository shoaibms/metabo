"""
Transformation Evaluation Script for Metabolomics Data Analysis
------------------------------------------------------------

This script evaluates the effectiveness of data transformations (specifically asinh)
by generating two key visualisations:

1. MA Plot (Mean-Difference Plot):
   - Visualises intensity-dependent ratio of two datasets.
   - X-axis: Average of log intensities (A).
   - Y-axis: Difference of log intensities (M).
   - Helps identify intensity-dependent bias.

2. Relative Standard Deviation (RSD) Plot:
   - Compares data variability before and after transformation.
   - Uses violin plots to show the full distribution.
   - Includes swarm plots for individual data points.

The analysis focuses on comparing original vs asinh-transformed data to assess
transformation effectiveness in improving data properties for metabolomics data.

Required Libraries:
- pandas: Data manipulation
- numpy: Numerical operations
- matplotlib: Basic plotting
- seaborn: Statistical visualisation
- statsmodels: LOWESS smoothing for MA plots
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from typing import Dict, Any, Tuple, List

# --- Configuration for Plots ---
# Using a dictionary for plot settings to make them easy to manage.
PLOT_STYLE_CONFIG = {
    'ma_plot': {
        'cmap': 'Greens',
        'title_size': 16,
        'label_size': 15,
        'tick_size': 14,
        'legend_size': 14,
        'cbar_label_size': 14,
    },
    'rsd_plot': {
        'before_color': '#f1facf',
        'after_color': '#68f7d8',
        'title_size': 18,
        'label_size': 17,
        'tick_size': 16,
    }
}


def load_data(file_path: str) -> pd.DataFrame:
    """
    Load data from a CSV file.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        pd.DataFrame: The loaded data as a pandas DataFrame.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file was not found at: {file_path}")
    return pd.read_csv(file_path)


def extract_n_cluster_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract columns containing 'N_Cluster' in their names.

    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: A subset of the DataFrame containing only N_Cluster columns.
    """
    return df.filter(regex="N_Cluster")


def _plot_single_ma(ax: plt.Axes, A: np.ndarray, M: np.ndarray, title: str, cmap: str,
                    style_config: Dict[str, Any], cbar_label: str):
    """Helper function to plot a single MA plot."""
    scatter = ax.scatter(A, M, alpha=0.5, c=M, cmap=cmap)
    
    # Add LOWESS trend line
    lowess = sm.nonparametric.lowess(M, A, frac=0.3)
    ax.plot(lowess[:, 0], lowess[:, 1], color='red', linewidth=2)
    
    ax.set_title(title, fontsize=style_config['title_size'])
    ax.set_xlabel('A (average log-intensity)', fontsize=style_config['label_size'])
    ax.set_ylabel('M (log ratio)', fontsize=style_config['label_size'])
    ax.tick_params(axis='both', which='major', labelsize=style_config['tick_size'])
    
    cbar = plt.gcf().colorbar(scatter, ax=ax)
    cbar.ax.tick_params(labelsize=style_config['legend_size'])
    cbar.set_label(cbar_label, fontsize=style_config['cbar_label_size'])


def plot_ma_transform(before_df: pd.DataFrame, after_df: pd.DataFrame, title: str):
    """
    Generate MA plots to compare data before and after transformation.

    MA plots (Bland-Altman plots) are used to visualise the intensity-dependent
    ratio of two datasets, which is helpful for identifying systematic bias.

    Args:
        before_df (pd.DataFrame): The original data.
        after_df (pd.DataFrame): The transformed data.
        title (str): The main title for the plot comparison.
    """
    style_config = PLOT_STYLE_CONFIG['ma_plot']
    
    # Add 1 to avoid log(0) issues
    before = before_df + 1
    after = after_df + 1

    # Ensure there are at least two columns for comparison
    if before.shape[1] < 2 or after.shape[1] < 2:
        raise ValueError("Input DataFrames must have at least two columns for MA plot.")

    # Calculate M (log ratio) and A (average log intensity)
    M_before = np.log2(before.iloc[:, 0]) - np.log2(before.iloc[:, 1])
    A_before = 0.5 * (np.log2(before.iloc[:, 0]) + np.log2(before.iloc[:, 1]))

    M_after = np.log2(after.iloc[:, 0]) - np.log2(after.iloc[:, 1])
    A_after = 0.5 * (np.log2(after.iloc[:, 0]) + np.log2(after.iloc[:, 1]))

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    _plot_single_ma(axes[0], A_before, M_before,
                    f'MA-plot Before Transformation - {title}',
                    style_config['cmap'], style_config, 'M Before')

    _plot_single_ma(axes[1], A_after, M_after,
                    f'MA-plot After Transformation - {title}',
                    style_config['cmap'], style_config, 'M After')

    plt.tight_layout()
    plt.show()


def calculate_rsd(data: pd.DataFrame) -> np.ndarray:
    """
    Calculate the Relative Standard Deviation (RSD) for each feature.

    RSD = (Standard Deviation / Mean) * 100

    Args:
        data (pd.DataFrame): The input data, with features as columns.

    Returns:
        np.ndarray: An array of RSD values.
    """
    mean = np.mean(data, axis=0)
    std_dev = np.std(data, axis=0)
    
    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        rsd = np.where(mean != 0, (std_dev / mean) * 100, np.nan)
    return rsd


def plot_rsd_violin(before_df: pd.DataFrame, after_df: pd.DataFrame, title: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate violin plots comparing RSD distributions before and after transformation.

    Args:
        before_df (pd.DataFrame): The original data.
        after_df (pd.DataFrame): The transformed data.
        title (str): The plot title.

    Returns:
        tuple: A tuple containing RSD values for 'before' and 'after' data.
    """
    style_config = PLOT_STYLE_CONFIG['rsd_plot']
    
    rsd_before = calculate_rsd(before_df)
    rsd_after = calculate_rsd(after_df)

    results_df = pd.DataFrame({
        'RSD': np.concatenate([rsd_before, rsd_after]),
        'Condition': ['Before'] * len(rsd_before) + ['After'] * len(rsd_after)
    })

    results_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    results_df.dropna(inplace=True)

    plt.figure(figsize=(8, 6))
    sns.violinplot(
        x='Condition', y='RSD', data=results_df,
        inner=None, palette=[style_config['before_color'], style_config['after_color']]
    )
    sns.swarmplot(x='Condition', y='RSD', data=results_df, color='k', alpha=0.4, size=4)

    plt.title(f'Relative Standard Deviation (RSD) - {title}', fontsize=style_config['title_size'])
    plt.xlabel('Data Type', fontsize=style_config['label_size'])
    plt.ylabel('RSD (%)', fontsize=style_config['label_size'])
    plt.xticks(fontsize=style_config['tick_size'])
    plt.yticks(fontsize=style_config['tick_size'])
    plt.tight_layout()
    plt.show()

    return rsd_before, rsd_after


def save_metrics_to_csv(metrics: Dict[str, List], file_path: str):
    """
    Save calculated metrics to a CSV file.

    Args:
        metrics (dict): A dictionary containing metric values.
        file_path (str): The path for the output CSV file.
    """
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(file_path, index=False)
    print(f"Metrics successfully saved to {file_path}")


def main():
    """Main function to execute the data transformation evaluation."""
    # As requested, file paths are hardcoded.
    # For broader use, consider making these configurable (e.g., command-line args).
    file_paths = {
        'original': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if.csv",
        'asinh': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_asinh.csv"
    }
    
    try:
        original_df = load_data(file_paths['original'])
        original_n_cluster = extract_n_cluster_data(original_df)

        transformed_df = load_data(file_paths['asinh'])
        transformed_n_cluster = extract_n_cluster_data(transformed_df)
        
        transformation_name = 'asinh'
        title = transformation_name.replace('_', ' ').title()

        # Generate plots
        plot_ma_transform(original_n_cluster, transformed_n_cluster, title)
        rsd_before, rsd_after = plot_rsd_violin(original_n_cluster, transformed_n_cluster, title)

        # Calculate and store metrics
        # Note: Corrected metric names to accurately reflect the calculation (mean RSD).
        # The original 'rMAD_Percentage_of_Median' seemed incorrect.
        metrics = {
            'Transformation': [title],
            'Mean_RSD_Before': [np.nanmean(rsd_before)],
            'Mean_RSD_After': [np.nanmean(rsd_after)],
        }
        
        # Define output path and save metrics
        output_dir = os.path.dirname(file_paths['original'])
        output_path = os.path.join(output_dir, 'transformation_metrics.csv')
        save_metrics_to_csv(metrics, output_path)

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except ValueError as e:
        print(f"Data Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()
