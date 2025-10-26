"""
Transformation Evaluation Script for Metabolomics Data
---------------------------------------------------

This script evaluates the effectiveness of data transformations in metabolomics analysis
using MA plots and Relative Standard Deviation (RSD) visualisation approaches.

Key Visualisations:
1. MA Plots: Assess intensity-dependent bias in the data
   - M (log ratio) vs A (average log-intensity)
   - LOWESS smoothing curve shows systematic trends
   
2. Violin Plots: Compare RSD distributions before/after transformation
   - Shows full distribution shape with density estimation
   - Includes swarm overlay for individual data points
   
Methods implemented:
- MA transformation for paired sample comparison
- RSD calculation for measurement precision assessment
- LOWESS smoothing for trend visualisation
- Violin plots with swarm overlay for distribution comparison

Metrics saved:
- Transformation method
- RSD before/after transformation
- rMAD percentage of median before/after
"""
from typing import Dict, Any, Tuple
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

FONT_SETTINGS = {
    'title': 16,
    'label': 14,
    'tick': 12
}

def load_data(file_path: Path) -> pd.DataFrame:
    """Load dataset from CSV file."""
    return pd.read_csv(file_path)

def extract_n_cluster_data(df: pd.DataFrame) -> pd.DataFrame:
    """Extract columns containing N_Cluster pattern."""
    return df.filter(regex="N_Cluster")

def plot_ma_transform(before_df: pd.DataFrame, after_df: pd.DataFrame, title: str, cmap: str = 'Greens') -> None:
    """
    Create MA plots comparing data before and after transformation.
    
    MA plots are used to visualize intensity-dependent ratio of changes:
    - M (y-axis): log ratio between conditions
    - A (x-axis): average log intensity
    
    Args:
        before_df: DataFrame before transformation
        after_df: DataFrame after transformation
        title: Plot title
        cmap: Colormap for scatter points
    """
    # Add 1 to avoid log(0)
    before = before_df + 1
    after = after_df + 1

    # Calculate M and A values
    M_before = np.log2(before.iloc[:, 0]) - np.log2(before.iloc[:, 1])
    A_before = 0.5 * (np.log2(before.iloc[:, 0]) + np.log2(before.iloc[:, 1]))
    M_after = np.log2(after.iloc[:, 0]) - np.log2(after.iloc[:, 1])
    A_after = 0.5 * (np.log2(after.iloc[:, 0]) + np.log2(after.iloc[:, 1]))

    # Create subplot figure
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))

    # Plot before transformation
    scatter1 = axes[0].scatter(A_before, M_before, alpha=0.5, c=M_before, cmap=cmap)
    lowess_before = sm.nonparametric.lowess(M_before, A_before, frac=0.3)
    axes[0].plot(lowess_before[:, 0], lowess_before[:, 1], color='red')
    axes[0].set_title(f'MA-transform Before Transformation - {title}', fontsize=FONT_SETTINGS['title'])
    axes[0].set_xlabel('A (average log-intensity)', fontsize=FONT_SETTINGS['label'])
    axes[0].set_ylabel('M (log ratio)', fontsize=FONT_SETTINGS['label'])
    axes[0].tick_params(axis='both', which='major', labelsize=FONT_SETTINGS['tick'])
    fig.colorbar(scatter1, ax=axes[0], label='M Before')

    # Plot after transformation
    scatter2 = axes[1].scatter(A_after, M_after, alpha=0.5, c=M_after, cmap=cmap)
    lowess_after = sm.nonparametric.lowess(M_after, A_after, frac=0.3)
    axes[1].plot(lowess_after[:, 0], lowess_after[:, 1], color='red')
    axes[1].set_title(f'MA-transform After Transformation - {title}', fontsize=FONT_SETTINGS['title'])
    axes[1].set_xlabel('A (average log-intensity)', fontsize=FONT_SETTINGS['label'])
    axes[1].set_ylabel('M (log ratio)', fontsize=FONT_SETTINGS['label'])
    axes[1].tick_params(axis='both', which='major', labelsize=FONT_SETTINGS['tick'])
    fig.colorbar(scatter2, ax=axes[1], label='M After')

    plt.tight_layout()
    plt.show()

def calculate_rsd(data: pd.DataFrame) -> pd.Series:
    """
    Calculate Relative Standard Deviation (RSD) for each column.
    
    RSD = (Standard Deviation / Mean) * 100
    """
    mean = np.mean(data, axis=0)
    std_dev = np.std(data, axis=0)
    
    # Handle division by zero
    rsd = np.divide(std_dev, mean, out=np.full_like(mean, np.nan), where=mean!=0) * 100
    return rsd

def calculate_rmad(series: pd.Series) -> float:
    """
    Calculate relative Median Absolute Deviation for a data series.
    
    rMAD = median(|x - median(x)|) / median(x)
    """
    median = np.median(series)
    if median == 0:
        return np.nan
    mad = np.median(np.abs(series - median))
    return mad / median

def plot_rsd_violin(before_df: pd.DataFrame, after_df: pd.DataFrame, title: str, before_color: str = '#649e68', 
                   after_color: str = '#b3e6b5') -> Tuple[pd.Series, pd.Series]:
    """
    Create violin plots comparing RSD distributions before and after transformation.
    
    Args:
        before_df, after_df: DataFrames before and after transformation
        title: Plot title
        before_color, after_color: Colors for before/after distributions
        
    Returns:
        Tuple of (rsd_before, rsd_after) Series
    """
    # Calculate RSD values
    rsd_before = calculate_rsd(before_df)
    rsd_after = calculate_rsd(after_df)

    # Prepare data for plotting
    results_df = pd.DataFrame({
        'RSD': np.concatenate([rsd_before, rsd_after]),
        'Condition': ['Before'] * len(rsd_before) + ['After'] * len(rsd_after)
    })

    # Remove infinite values and NaNs
    results_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    results_df.dropna(inplace=True)

    # Create violin plot
    plt.figure(figsize=(10, 6))
    sns.violinplot(x='Condition', y='RSD', data=results_df, inner=None, 
                  palette=[before_color, after_color])
    sns.swarmplot(x='Condition', y='RSD', data=results_df, color='k', alpha=0.5)
    
    # Set plot attributes
    plt.title(f'Violin Plot of Relative Standard Deviation (RSD) - {title}', 
             fontsize=FONT_SETTINGS['title'])
    plt.xlabel('Condition', fontsize=FONT_SETTINGS['label'])
    plt.ylabel('RSD (%)', fontsize=FONT_SETTINGS['label'])
    plt.xticks(fontsize=FONT_SETTINGS['tick'])
    plt.yticks(fontsize=FONT_SETTINGS['tick'])
    plt.show()

    return rsd_before, rsd_after

def save_metrics_to_csv(metrics: Dict[str, Any], file_path: Path) -> None:
    """Save transformation metrics to CSV file."""
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(file_path, index=False)
    print(f"Metrics saved to {file_path}")

def main():
    """Main function to run the transformation evaluation."""
    # Define file paths for all transformations
    base_path = Path(r"C:\Users\ms\Desktop\data_chem\data\transformation")
    file_paths = {
        'original': base_path / "n_l_if.csv",
        'anscombe': base_path / "n_l_if_anscombe.csv",
        'asinh': base_path / "n_l_if_asinh.csv",
        'boxcox': base_path / "n_l_if_boxcox.csv",
        'glog': base_path / "n_l_if_glog.csv",
        'log': base_path / "n_l_if_log.csv",
        'sqrt': base_path / "n_l_if_sqrt.csv",
        'yeojohnson': base_path / "n_l_if_yeojohnson.csv"
    }

    # Load original data
    original_df = load_data(file_paths['original'])
    original_n_cluster = extract_n_cluster_data(original_df)
    
    # Initialise metrics dictionary
    metrics = {
        'Transformation': [],
        'RSD_Before': [],
        'RSD_After': [],
        'rMAD_Before': [],
        'rMAD_After': []
    }

    # Process each transformation
    for key, path in file_paths.items():
        if key == 'original':
            continue
        
        transformed_df = load_data(path)
        transformed_n_cluster = extract_n_cluster_data(transformed_df)
        title = key.replace('_', ' ').title()

        # Generate plots and calculate metrics
        plot_ma_transform(original_n_cluster, transformed_n_cluster, title)
        rsd_before, rsd_after = plot_rsd_violin(original_n_cluster, transformed_n_cluster, title)
        
        rmad_before = original_n_cluster.apply(calculate_rmad).mean()
        rmad_after = transformed_n_cluster.apply(calculate_rmad).mean()

        # Store metrics
        metrics['Transformation'].append(title)
        metrics['RSD_Before'].append(rsd_before.mean())
        metrics['RSD_After'].append(rsd_after.mean())
        metrics['rMAD_Before'].append(rmad_before)
        metrics['rMAD_After'].append(rmad_after)

    # Save results
    save_metrics_to_csv(metrics, base_path / 'transformation_metrics_2.csv')

if __name__ == "__main__":
    main()
