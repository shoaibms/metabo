"""
Logistic Regression Results Analysis and Visualization
---------------------------------------------------

This script analyzes logistic regression results from multiple metabolomics datasets
to identify significant associations between metabolites and experimental factors.

The analysis:
1. Processes results from four complementary datasets (n_l, n_r, p_l, p_r)
2. Calculates percentage of significant associations for each predictor
3. Generates visualizations comparing significance patterns across datasets
4. Creates both tabular and graphical summaries of the findings

Input files:
- CSV files containing logistic regression results for each dataset
- Expected columns for p-values: Genotype_G2, TMT_1, Day_2, Rep

Output:
- Summary CSV with percentages of significant associations
- Grouped bar plot visualization of the results
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Tuple, Any

def get_file_paths() -> Tuple[Dict[str, str], Dict[str, str]]:
    """Defines and returns paths for input and output files."""
    base_path = r'C:\Users\USER\Desktop\data_chem\data\data\LogisticRegression'
    
    input_files = {
        f'{base_path}\\n_column_data_l_clean_results.csv': 'n_l',
        f'{base_path}\\n_column_data_r_clean_results.csv': 'n_r',
        f'{base_path}\\p_column_data_l_clean_results.csv': 'p_l',
        f'{base_path}\\p_column_data_r_clean_results.csv': 'p_r'
    }
    
    output_files = {
        'summary_table': f'{base_path}\\summary_table_significant_associations.csv',
        'grouped_bar_plot': f'{base_path}\\grouped_bar_significant_associations_datasets.png'
    }
    
    return input_files, output_files

def calculate_significant_associations(df: pd.DataFrame, column: str, threshold: float = 0.05) -> int:
    """
    Calculates the number of significant associations for a predictor column.

    A value is considered significant if it is less than the threshold.
    Non-numeric values are ignored.

    Args:
        df: DataFrame containing regression results.
        column: Name of the predictor column (containing p-values).
        threshold: Significance threshold.

    Returns:
        The number of significant associations.
    """
    if column not in df.columns:
        return 0
    
    p_values = pd.to_numeric(df[column], errors='coerce')
    return (p_values < threshold).sum()

def process_dataset(file_path: str, dataset_name: str) -> Dict[str, Any] | None:
    """
    Processes a single dataset to calculate summary statistics.

    Args:
        file_path: Path to the dataset file.
        dataset_name: Name identifier for the dataset.

    Returns:
        A dictionary with summary statistics, or None if the file is not found.
    """
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Warning: File not found at {file_path}. Skipping.")
        return None

    df_length = len(df)
    if df_length == 0:
        return {'Dataset': dataset_name, 'Genotype': 0, 'Treatment': 0, 'Day': 0, 'Replicate': 0}

    return {
        'Dataset': dataset_name,
        'Genotype': (calculate_significant_associations(df, 'Genotype_G2') / df_length) * 100,
        'Treatment': (calculate_significant_associations(df, 'TMT_1') / df_length) * 100,
        'Day': (calculate_significant_associations(df, 'Day_2') / df_length) * 100,
        'Replicate': (calculate_significant_associations(df, 'Rep') / df_length) * 100
    }

def create_summary_visualization(summary_df: pd.DataFrame, output_path: str, font_size: int = 21):
    """
    Creates and saves a grouped bar plot visualization from summary data.

    Args:
        summary_df: DataFrame with summary statistics.
        output_path: Path to save the visualization.
        font_size: Base font size for plot elements.
    """
    summary_melted = summary_df.melt(
        id_vars='Dataset',
        var_name='Predictor',
        value_name='Percentage'
    )
    
    color_palette = sns.color_palette("Greens", len(summary_melted['Predictor'].unique()))
    
    plt.figure(figsize=(8, 5.85))
    sns.barplot(
        x='Dataset',
        y='Percentage',
        hue='Predictor',
        data=summary_melted,
        palette=color_palette
    )
    
    plt.xlabel('Dataset', fontsize=font_size)
    plt.ylabel('% of Significant Associations', fontsize=font_size)
    plt.xticks(fontsize=font_size - 1)
    plt.yticks(fontsize=font_size - 1)
    plt.legend(
        title='Predictor',
        fontsize=font_size - 2,
        title_fontsize=font_size,
        framealpha=0.6
    )
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

def main():
    """Main function to run the analysis and visualization workflow."""
    input_files, output_files = get_file_paths()
    
    summary_data = []
    for file_path, dataset_name in input_files.items():
        summary = process_dataset(file_path, dataset_name)
        if summary:
            summary_data.append(summary)
    
    if not summary_data:
        print("No data processed. Exiting.")
        return

    summary_df = pd.DataFrame(summary_data)
    
    summary_file_path = output_files['summary_table']
    summary_df.to_csv(summary_file_path, index=False)
    print(f"Summary table saved to: {summary_file_path}")
    
    plot_path = output_files['grouped_bar_plot']
    create_summary_visualization(summary_df, plot_path)
    print(f"Grouped bar chart saved to: {plot_path}")

if __name__ == "__main__":
    main()
