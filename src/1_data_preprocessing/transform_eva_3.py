"""
Data Transformation Evaluation for Metabolomics Analysis
======================================================

This script provides a framework for evaluating the impact of data transformations
on metabolomics datasets. It generates facet grid visualizations to compare the
distribution of two key statistical metrics—Coefficient of Variation (CV) and
Relative Median Absolute Deviation (rMAD)—before and after a specified
transformation.

The primary output consists of violin plots that allow for a direct visual
assessment of how a transformation affects data variability and dispersion.

Inputs:
-------
- An original data file in CSV format.
- A transformed data file in CSV format.
- Both files are expected to contain columns with the 'N_Cluster' identifier,
  which are used for the analysis.

Outputs:
--------
- A set of facet grid plots displaying:
  1. A comparison of CV distributions for original vs. transformed data.
  2. A comparison of rMAD distributions for original vs. transformed data.

"""

from pathlib import Path
from typing import Tuple, Dict, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Define font sizes for plots as a global constant
FONT_SETTINGS = {
    'title': 18,
    'xlabel': 16,
    'ylabel': 16,
    'ticks': 14,
    'legend': 14
}

def extract_n_cluster_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts columns containing 'N_Cluster' in their names.

    Args:
        df: Input DataFrame.

    Returns:
        A DataFrame containing only the 'N_Cluster' columns.
    """
    return df.filter(regex="N_Cluster")


def calculate_cv(series: pd.Series) -> float:
    """
    Calculates the Coefficient of Variation (CV) for a data series.

    Formula:
        CV = (standard deviation / mean)

    Args:
        series: Input data series.

    Returns:
        The CV value. Returns NaN if the mean is zero.
    """
    mean = np.mean(series)
    if mean == 0:
        return np.nan
    return np.std(series) / mean


def calculate_rmad(series: pd.Series) -> float:
    """
    Calculates the Relative Median Absolute Deviation (rMAD) for a data series.

    Formula:
        rMAD = median(|x - median(x)|) / median(x)

    Args:
        series: Input data series.

    Returns:
        The rMAD value. Returns NaN if the median is zero.
    """
    median = np.median(series)
    if median == 0:
        return np.nan
    mad = np.median(np.abs(series - median))
    return mad / median


def calculate_metrics(df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    """
    Calculates CV and rMAD for all columns in a DataFrame.

    Args:
        df: Input DataFrame with features as columns.

    Returns:
        A tuple containing the CV and rMAD values for each column.
    """
    cv_values = df.apply(calculate_cv)
    rmad_values = df.apply(calculate_rmad)
    return cv_values, rmad_values


def load_data(original_path: Path, transformed_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads the original and transformed datasets from CSV files.

    Args:
        original_path: Path to the original data CSV file.
        transformed_path: Path to the transformed data CSV file.

    Returns:
        A tuple containing the original and transformed DataFrames.
        
    Raises:
        FileNotFoundError: If a file is not found at the specified path.
        pd.errors.EmptyDataError: If the CSV file is empty.
    """
    try:
        original_df = pd.read_csv(original_path)
        transformed_df = pd.read_csv(transformed_path)
        return original_df, transformed_df
    except (FileNotFoundError, pd.errors.EmptyDataError) as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        raise

def prepare_data_for_plotting(
    original_df: pd.DataFrame, 
    transformed_df: pd.DataFrame,
    transformation_name: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Prepares data for visualization by calculating metrics and formatting.

    Args:
        original_df: The original data.
        transformed_df: The transformed data.
        transformation_name: The name of the transformation applied (e.g., 'asinh').

    Returns:
        A tuple of DataFrames for CV and rMAD, ready for plotting.
    """
    orig_n_cluster = extract_n_cluster_data(original_df)
    trans_n_cluster = extract_n_cluster_data(transformed_df)

    orig_cv, orig_rmad = calculate_metrics(orig_n_cluster)
    trans_cv, trans_rmad = calculate_metrics(trans_n_cluster)

    def _prepare_metric_df(orig_metric: pd.Series, trans_metric: pd.Series, metric_name: str) -> pd.DataFrame:
        """Helper to format metric data for plotting."""
        df = pd.DataFrame({'Original': orig_metric, 'Transformed': trans_metric})
        df = df.reset_index().melt(id_vars='index', var_name='Type', value_name=metric_name)
        df['Transformation'] = transformation_name
        return df

    cv_df = _prepare_metric_df(orig_cv, trans_cv, 'CV')
    rmad_df = _prepare_metric_df(orig_rmad, trans_rmad, 'rMAD')
    
    return cv_df, rmad_df


def plot_facet_grid(df: pd.DataFrame, metric: str, font_settings: Dict[str, Any]):
    """
    Creates a facet grid violin plot to compare distributions.

    Args:
        df: Data prepared for plotting.
        metric: The name of the metric to plot (e.g., 'CV' or 'rMAD').
        font_settings: A dictionary of font sizes for plot elements.
    """
    g = sns.FacetGrid(df, col="Transformation", sharex=False, sharey=False, height=4)
    g.map_dataframe(sns.violinplot, x='Type', y=metric, palette='BuGn')
    
    g.set_axis_labels("Data Type", metric)
    g.set_titles(col_template="{col_name}", size=font_settings['title'])

    for ax in g.axes.flat:
        if ax:
            ax.set_xlabel('Data Type', fontsize=font_settings['xlabel'])
            ax.set_ylabel(metric, fontsize=font_settings['ylabel'])
            ax.tick_params(axis='both', which='major', labelsize=font_settings['ticks'])

    g.fig.subplots_adjust(top=0.85, hspace=0.4)
    g.fig.suptitle(f'Violin Plot of {metric}', fontsize=font_settings['title'])
    
    plt.tight_layout()
    plt.show()


def main():
    """
    Main function to run the data transformation evaluation.
    """
    original_path = Path(r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if.csv")
    transformed_path = Path(r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_asinh.csv")
    transformation_name = 'asinh'
    
    try:
        original_df, transformed_df = load_data(original_path, transformed_path)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return

    cv_df, rmad_df = prepare_data_for_plotting(
        original_df, transformed_df, transformation_name
    )
    
    plot_facet_grid(cv_df, 'CV', FONT_SETTINGS)
    plot_facet_grid(rmad_df, 'rMAD', FONT_SETTINGS)


if __name__ == "__main__":
    main()
