"""
Relative Median Absolute Deviation (rMAD) Based Variable Selection
---------------------------------------------------------------

This script performs variable selection in metabolomics data based on the relative Median 
Absolute Deviation (rMAD) criterion. Variables with high rMAD values, indicating high 
variability relative to their median, are removed from the dataset.

Key features:
- Calculates rMAD for metabolite clusters (N_Cluster and P_Cluster variables)
- Removes variables above a specified rMAD threshold
- Preserves non-cluster variables
- Generates both cleaned datasets and lists of removed variables

Input:
- CSV files containing metabolomics data with N_Cluster or P_Cluster columns
- rMAD threshold for variable removal (default: 30%)

Output:
- Cleaned CSV files with high-rMAD variables removed
- CSV files listing removed variables for documentation
"""

import os
from typing import List, Tuple

import numpy as np
import pandas as pd

# --- Configuration ---
# RMAD threshold for variable removal. Features with rMAD greater than this will be removed.
RMAD_THRESHOLD = 30.0


def calculate_rmad(data: pd.DataFrame) -> pd.Series:
    """
    Calculate Relative Median Absolute Deviation (rMAD) for each column.
    
    rMAD = (median(|x - median(x)|) / median(x)) * 100
    
    Args:
        data: DataFrame containing numeric columns.
        
    Returns:
        Series containing rMAD values for each column.
    """
    median = data.median()
    # Handle potential division by zero if median is 0
    median_no_zero = median.replace(0, np.nan)
    mad = (np.abs(data - median)).median()
    rmad = (mad / median_no_zero) * 100
    return rmad.fillna(0)


def _get_column_types(data: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """
    Identify and separate metabolite cluster columns from other columns.

    Args:
        data: Input DataFrame.

    Returns:
        A tuple containing two lists:
        - List of column names starting with 'N_Cluster' or 'P_Cluster'.
        - List of all other column names.
    """
    cluster_cols = [
        col for col in data.columns if col.startswith(("N_Cluster", "P_Cluster"))
    ]
    non_cluster_cols = [col for col in data.columns if col not in cluster_cols]
    return cluster_cols, non_cluster_cols


def filter_variables_by_rmad(
    data: pd.DataFrame,
    rmad_values: pd.Series,
    non_cluster_cols: List[str],
    threshold: float,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Filter variables based on the rMAD threshold, keeping non-cluster columns.

    Args:
        data: The original input DataFrame.
        rmad_values: Series containing rMAD values for cluster columns.
        non_cluster_cols: A list of columns to always keep.
        threshold: The maximum acceptable rMAD value.

    Returns:
        A tuple containing:
        - DataFrame with high-rMAD variables removed.
        - List of removed column names.
    """
    columns_to_keep = rmad_values[rmad_values <= threshold].index.tolist()
    columns_to_remove = rmad_values[rmad_values > threshold].index.tolist()

    final_columns = non_cluster_cols + columns_to_keep
    clean_data = data[final_columns].copy()

    return clean_data, columns_to_remove


def save_results(
    clean_data: pd.DataFrame, removed_vars: List[str], original_path: str
) -> Tuple[str, str]:
    """
    Save the cleaned data and the list of removed variables to new files.

    Args:
        clean_data: The filtered DataFrame.
        removed_vars: A list of removed column names.
        original_path: The path to the original input file.

    Returns:
        A tuple containing the paths to the two saved files.
    """
    directory = os.path.dirname(original_path)
    filename = os.path.basename(original_path)

    clean_data_path = os.path.join(directory, f"cleaned_{filename}")
    clean_data.to_csv(clean_data_path, index=False)

    removed_vars_path = os.path.join(directory, f"removed_variables_{filename}")
    pd.DataFrame(removed_vars, columns=["Removed_Variables"]).to_csv(
        removed_vars_path, index=False
    )

    return clean_data_path, removed_vars_path


def process_file(file_path: str, rmad_threshold: float) -> None:
    """
    Process a single file to remove variables with high rMAD values.

    Args:
        file_path: The path to the input CSV file.
        rmad_threshold: The maximum acceptable rMAD value.
    """
    try:
        data = pd.read_csv(file_path)

        cluster_columns, non_cluster_columns = _get_column_types(data)
        cluster_data = data[cluster_columns]

        rmad_values = calculate_rmad(cluster_data)
        clean_data, removed_variables = filter_variables_by_rmad(
            data, rmad_values, non_cluster_columns, rmad_threshold
        )

        clean_path, removed_path = save_results(
            clean_data, removed_variables, file_path
        )

        print(f"\nResults for {os.path.basename(file_path)}:")
        print(f"  - Cleaned data saved to: {clean_path}")
        print(f"  - List of removed variables saved to: {removed_path}")
        print(f"  - Number of variables removed: {len(removed_variables)}")

    except FileNotFoundError:
        print(f"Error: The file was not found at {file_path}")
    except Exception as e:
        print(f"An error occurred while processing {file_path}: {e}")


def main():
    """
    Main function to execute the rMAD-based variable selection process.
    """
    # NOTE: Please update the list below with the full paths to your input files.
    file_paths = [
        r"C:\Users\ms\Desktop\data_chem\data\CV_rMAD\n_l_if_asinh.csv",
        r"C:\Users\ms\Desktop\data_chem\data\CV_rMAD\n_r_if_asinh.csv",
        r"C:\Users\ms\Desktop\data_chem\data\CV_rMAD\p_l_if_asinh.csv",
        r"C:\Users\ms\Desktop\data_chem\data\CV_rMAD\p_r_if_asinh.csv",
    ]

    for file_path in file_paths:
        if os.path.exists(file_path):
            process_file(file_path, RMAD_THRESHOLD)
        else:
            print(f"Warning: The file path does not exist, skipping: {file_path}")


if __name__ == "__main__":
    main()
