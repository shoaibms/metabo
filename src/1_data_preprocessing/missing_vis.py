"""
Visualisation of Missing Data in Metabolomics Data
==================================================

This script provides a visualisation of missing values in metabolomics data,
focusing on variables prefixed with 'N_Cluster_'. It generates a heatmap
to represent the data's completeness.

In the generated heatmap:
- Missing values are colored green.
- Non-missing (present) values are colored light gray.
- The legend displays the percentage of missing and non-missing values.

This visualisation is useful for:
- Understanding the distribution and pattern of missing data.
- Identifying if the missingness is systematic or random.
- Quantifying the extent of missing data across features and samples.
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from pathlib import Path
from typing import Dict, Any

# --- Configuration ---
CONFIG = {
    'data_path': Path(r'C:\Users\ms\Desktop\data_chem\data\old_2\n_column_data_l_deducted.csv'),
    'output_file': Path('figure/missing_data_visualization.png'),
    'font_size': 15,
    'cluster_prefix': 'N_Cluster_',
    'colors': {
        'non_missing': '#e1e6eb',  # Light gray
        'missing': '#2ca02c'       # Green
    }
}

def plot_missing_data_matrix(df: pd.DataFrame, config: Dict[str, Any]) -> None:
    """
    Creates and saves a heatmap visualization of missing data patterns.

    Args:
        df (pd.DataFrame): Input DataFrame containing metabolomics data.
        config (Dict[str, Any]): A dictionary containing configuration parameters.
    """
    plt.rcParams.update({'font.size': config['font_size']})

    # Filter for columns of interest
    n_cluster_data = df.loc[:, df.columns.str.startswith(config['cluster_prefix'])]

    # Create a boolean mask for missing values (True where value is missing)
    missing_mask = n_cluster_data.isnull()

    # Calculate statistics on missing data
    total_elements = missing_mask.size
    missing_elements = missing_mask.sum().sum()
    
    if total_elements == 0:
        print("No data to plot.")
        return
        
    missing_percentage = (missing_elements / total_elements) * 100
    non_missing_percentage = 100 - missing_percentage

    # Define colors for visualization
    custom_cmap = ListedColormap([config['colors']['non_missing'], config['colors']['missing']])

    # Create and configure the plot
    fig, ax = plt.subplots(figsize=(5, 4))
    
    # Remove plot borders for a cleaner appearance
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Display the missing data mask as an image
    ax.imshow(missing_mask, 
              aspect='auto',
              interpolation='none',
              cmap=custom_cmap)

    # Configure axes with dynamically generated labels
    num_entries, num_variables = n_cluster_data.shape
    ax.set_xlabel(f'Variables (N={num_variables})')
    ax.set_ylabel(f'Entry (N={num_entries})')
    
    # Remove tick marks to reduce visual clutter
    ax.set_xticks([])
    ax.set_yticks([])

    # Create and configure the legend
    legend_labels = [
        f'Non-Missing ({non_missing_percentage:.2f}%)',
        f'Missing ({missing_percentage:.2f}%)'
    ]
    patches = [Patch(color=config['colors'][key], label=label)
               for key, label in zip(['non_missing', 'missing'], legend_labels)]
    
    leg = ax.legend(handles=patches,
                    loc='upper right',
                    frameon=True)
    leg.get_frame().set_edgecolor('none')

    plt.tight_layout()
    
    # Save the figure
    output_path = config['output_file']
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")
    plt.close(fig)

def main(config: Dict[str, Any]) -> None:
    """
    Main function to load data and generate the missing data visualization.
    """
    try:
        data = pd.read_csv(config['data_path'])
    except FileNotFoundError:
        print(f"Error: The file was not found at {config['data_path']}")
        return

    # Generate and save the plot
    plot_missing_data_matrix(data, config)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize missing data in metabolomics data.")
    parser.add_argument("--data_path", type=Path, default=CONFIG['data_path'],
                        help="Path to the input CSV file.")
    parser.add_argument("--output_file", type=Path, default=CONFIG['output_file'],
                        help="Path to save the output plot.")
    args = parser.parse_args()

    current_config = CONFIG.copy()
    current_config['data_path'] = args.data_path
    current_config['output_file'] = args.output_file
    
    main(current_config)
