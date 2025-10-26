"""
Imputation Validation Visualization Script
----------------------------------------

This script generates pairwise box plots to visually compare an original metabolomics 
dataset with its imputed counterpart. It focuses on a random subset of features 
(variables containing 'N_Cluster' in their names) to assess how well the 
imputation method preserves the original data distribution.

The script performs the following steps:
1.  Loads the original (with missing values) and imputed datasets.
2.  Randomly selects a specified number of 'N_Cluster' columns for comparison.
3.  For visualization purposes, it fills any missing values in the selected columns of the 
    original data using the median of each respective column.
4.  Restructures both original and imputed data into a long format suitable for plotting.
5.  Generates and displays a publication-quality box plot comparing the distributions.
6.  Saves the generated plot to a file.

This process aids in validating the imputation by providing a clear visual
representation of the data's characteristics before and after imputation.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict, Any

def load_datasets(original_path: str, imputed_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load original and imputed datasets from CSV files.
    
    Args:
        original_path (str): The file path to the original dataset.
        imputed_path (str): The file path to the imputed dataset.
    
    Returns:
        A tuple containing the original and imputed pandas DataFrames.
    """
    original_data = pd.read_csv(original_path)
    imputed_data = pd.read_csv(imputed_path)
    return original_data, imputed_data

def select_random_clusters(data: pd.DataFrame, n_clusters: int, random_seed: int) -> List[str]:
    """
    Randomly select a specified number of 'N_Cluster' columns from the dataset.
    
    Args:
        data (pd.DataFrame): The input DataFrame.
        n_clusters (int): The number of cluster columns to select.
        random_seed (int): The seed for the random number generator for reproducibility.
    
    Returns:
        A list of the selected column names.
    """
    np.random.seed(random_seed)
    cluster_columns = [col for col in data.columns if 'N_Cluster' in col]
    if len(cluster_columns) < n_clusters:
        raise ValueError(f"Not enough 'N_Cluster' columns to select {n_clusters}. Found {len(cluster_columns)}.")
    return np.random.choice(cluster_columns, size=n_clusters, replace=False).tolist()

def prepare_data_for_plotting(original_data: pd.DataFrame, imputed_data: pd.DataFrame, selected_columns: List[str]) -> pd.DataFrame:
    """
    Prepares the data for plotting by selecting, transforming, and combining subsets.

    For the original data subset, missing values are filled with the column medians
    to allow for visualization. Both subsets are then melted and combined.
    
    Args:
        original_data (pd.DataFrame): The original dataset.
        imputed_data (pd.DataFrame): The imputed dataset.
        selected_columns (list): A list of column names to be processed.
    
    Returns:
        A pandas DataFrame ready for plotting, containing combined data from both
        original and imputed sources.
    """
    # Extract and process original data
    original_subset = original_data[selected_columns].copy()
    original_subset.fillna(original_subset.median(), inplace=True)
    original_melted = original_subset.melt(var_name='Variable', value_name='Value')
    original_melted['Type'] = 'Original'
    
    # Process imputed data
    imputed_subset = imputed_data[selected_columns].copy()
    imputed_melted = imputed_subset.melt(var_name='Variable', value_name='Value')
    imputed_melted['Type'] = 'Imputed'
    
    # Combine datasets
    return pd.concat([original_melted, imputed_melted])

def create_comparison_boxplot(data: pd.DataFrame, palette: Dict[str, str], plot_params: Dict[str, Any], save_path: Path) -> None:
    """
    Creates, customizes, and saves a boxplot comparing original and imputed values.
    
    Args:
        data (pd.DataFrame): The prepared data for plotting.
        palette (dict): The color palette for different data types.
        plot_params (dict): Parameters for customizing the plot's appearance.
        save_path (Path): The file path to save the plot.
    """
    plt.figure(figsize=plot_params['figsize'])
    
    sns.boxplot(
        data=data,
        x='Variable',
        y='Value',
        hue='Type',
        palette=palette,
        boxprops=dict(edgecolor='none')
    )
    
    plt.xticks(rotation=90, fontsize=plot_params['tick_fontsize'])
    plt.yticks(fontsize=plot_params['tick_fontsize'])
    plt.xlabel('Variable', fontsize=plot_params['label_fontsize'])
    plt.ylabel('Value', fontsize=plot_params['label_fontsize'])
    plt.title(plot_params['title'], fontsize=plot_params['title_fontsize'])
    
    plt.legend(
        title='Data Type',
        title_fontsize=plot_params['legend_title_fontsize'],
        fontsize=plot_params['legend_fontsize']
    )
    
    plt.tight_layout()
    
    # Save the figure
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {save_path}")

    plt.show()

def main():
    """
    Main function to run the imputation validation and visualization pipeline.
    """
    # --- Configuration ---
    # File paths (kept as per user request)
    original_path = r"C:\Users\ms\Desktop\data_chem\data\old_2\n_column_data_l_deducted.csv"
    imputed_path = r"C:\Users\ms\Desktop\data_chem\imputated\complete imputed files\n_column_data_l_deducted_imputed_rf2.csv"
    output_dir = Path("figure")
    output_filename = "imputation_validation_boxplot.png"
    
    # Parameters
    N_CLUSTERS = 20
    RANDOM_SEED = 42
    
    # Color palette
    palette = {
        "Original": "#F0CACA",
        "Imputed": "#F3E2C4"
    }
    
    # Plot parameters
    plot_params = {
        'figsize': (20, 10),
        'tick_fontsize': 12,
        'label_fontsize': 14,
        'title_fontsize': 16,
        'legend_fontsize': 12,
        'legend_title_fontsize': 13,
        'title': f'Pairwise Box Plot of {N_CLUSTERS} Random N_Cluster Variables from Original and Imputed Data'
    }
    
    # --- Execution ---
    # Load and process data
    original_data, imputed_data = load_datasets(original_path, imputed_path)
    selected_columns = select_random_clusters(original_data, n_clusters=N_CLUSTERS, random_seed=RANDOM_SEED)
    combined_data = prepare_data_for_plotting(original_data, imputed_data, selected_columns)
    
    # Create and save visualization
    save_path = output_dir / output_filename
    create_comparison_boxplot(combined_data, palette, plot_params, save_path)

if __name__ == "__main__":
    main()
