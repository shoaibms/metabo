"""
Normality Test Visualization Script for Metabolomics Data Transformations
----------------------------------------------------------------------

This script creates density plots to visualize the results of normality tests
(Shapiro-Wilk and Anderson-Darling) across different data transformation methods.
It helps assess which transformation methods most effectively normalize the data,
saving the output plots as image files.

Visualizations:
1. Shapiro-Wilk test statistic distribution
2. Anderson-Darling test statistic distribution
3. Shapiro-Wilk p-value distribution
4. Anderson-Darling p-value distribution

Input:
- CSV files containing normality test results for different transformations.
- Each file should be named following a convention where the transformation
  method is identifiable, e.g., '..._{transformation_name}_normality_results.csv'.

Output:
- Four density plots saved as PNG files, comparing distributions across
  transformation methods.
"""
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# --- Configuration ---
CONFIG = {
    'input_dir': Path(r'C:\Users\USER\Desktop\data_chem\data\transformation'),
    'output_dir': Path('figure/normality_plots'),
    'hex_colors': [
        "#006400", "#228B22", "#32CD32", "#66CDAA", "#8FBC8F", "#98FB98",
        "#ADFF2F"
    ],
    'font_sizes': {
        'xlabel': 16,
        'ylabel': 16,
        'ticks': 14,
        'legend_title': 16,
        'legend_text': 15
    },
    'metrics_to_plot': [
        'Shapiro_Statistic', 'Anderson_Statistic', 'Shapiro_p_value',
        'Anderson_p_value'
    ]
}

def get_transformation_data(
        input_dir: Path) -> Tuple[pd.DataFrame, List[str]]:
    """
    Load and combine normality test results from multiple transformation files.

    Args:
        input_dir: A path to the directory containing CSV files with normality test results.

    Returns:
        A tuple containing:
        - A combined DataFrame with an added 'Transformation' column.
        - A list of transformation method names extracted from filenames.
    """
    dfs = []
    transformation_methods = []
    
    file_paths = list(input_dir.glob('*_normality_results.csv'))

    for file_path in file_paths:
        df = pd.read_csv(file_path)
        # Assumes filename format like: ..._{transformation_name}_normality_results.csv
        transformation_method = file_path.stem.split('_')[-3]

        df['Transformation'] = transformation_method
        dfs.append(df)
        transformation_methods.append(transformation_method)

    return pd.concat(dfs, ignore_index=True), sorted(list(set(transformation_methods)))


def create_density_plot(data: pd.DataFrame, x_column: str,
                        transformation_methods: List[str],
                        hex_colors: List[str], font_sizes: Dict[str, int],
                        output_path: Path):
    """
    Create and save a density plot for a specific normality test metric.

    Args:
        data: Combined normality test results.
        x_column: Column name to plot (e.g., 'Shapiro_Statistic').
        transformation_methods: Names of transformation methods to control plot order.
        hex_colors: Color palette for different transformations.
        font_sizes: Font size specifications.
        output_path: Path to save the generated plot.
    """
    plt.figure(figsize=(6, 4))

    sns.kdeplot(data=data,
                x=x_column,
                hue='Transformation',
                hue_order=transformation_methods,
                fill=True,
                common_norm=False,
                palette=hex_colors)

    plt.xlabel(x_column.replace('_', ' ').title(), fontsize=font_sizes['xlabel'])
    plt.ylabel('Density', fontsize=font_sizes['ylabel'])
    plt.xticks(fontsize=font_sizes['ticks'])
    plt.yticks(fontsize=font_sizes['ticks'])

    legend = plt.legend(
        title='Transformation Method',
        title_fontsize=font_sizes['legend_title'],
        fontsize=font_sizes['legend_text'],
    )
    legend.get_frame().set_alpha(0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved to {output_path}")


def main(config: Dict):
    """
    Main function to generate and save normality test visualization plots.
    """
    config['output_dir'].mkdir(parents=True, exist_ok=True)

    all_data, transformation_methods = get_transformation_data(config['input_dir'])

    for metric in config['metrics_to_plot']:
        output_path = config['output_dir'] / f"{metric}_density_plot.png"
        create_density_plot(all_data, metric, transformation_methods,
                            config['hex_colors'], config['font_sizes'], output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate density plots for normality test results.")
    parser.add_argument("--input_dir", type=Path, default=CONFIG['input_dir'],
                        help="Directory containing normality test CSV files.")
    parser.add_argument("--output_dir", type=Path, default=CONFIG['output_dir'],
                        help="Directory to save the output plots.")
    args = parser.parse_args()

    current_config = CONFIG.copy()
    current_config['input_dir'] = args.input_dir
    current_config['output_dir'] = args.output_dir
    
    main(current_config)
