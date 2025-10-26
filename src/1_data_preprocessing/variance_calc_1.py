"""
Metabolomics Data Quality Assessment via Variability Metrics.

This script evaluates the quality of metabolomics datasets by calculating two
key variability metrics for each feature:
1.  **Coefficient of Variation (CV)**: The ratio of the standard deviation to the
    mean, expressed as a percentage. It is a standardized measure of dispersion.
2.  **Relative Median Absolute Deviation (rMAD)**: A robust alternative to CV,
    calculated from the median absolute deviation and the median. It is less
    sensitive to outliers.

High values for these metrics (e.g., > 30%) can indicate low-quality features
that may need to be filtered or transformed prior to downstream analysis.

This script processes a predefined list of input CSV files, identifies columns
representing metabolite features (those prefixed with 'N_Cluster' or 'P_Cluster'),
calculates CV and rMAD, and reports the percentage of features exceeding
predefined thresholds for each file.

Usage:
    - Update the `FILE_PATHS` list with the paths to your data files.
    - Run the script from the command line: `python variance_calc_1.py`
"""

import pandas as pd
import numpy as np
import os
from typing import Tuple, List, Optional

# --- Configuration ---
# Quality thresholds (in percentage)
CV_THRESHOLD = 30.0
RMAD_THRESHOLD = 30.0

# NOTE: The file paths below are hardcoded as per user request.
# For general use, you may want to modify this to accept paths from a
# configuration file or command-line arguments.
FILE_PATHS = [
    r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\n_l_if_asinh.csv',
    r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\n_r_if_asinh.csv',
    r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\p_l_if_asinh.csv',
    r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\p_r_if_asinh.csv'
]
# --- End of Configuration ---

def calculate_cv(data: pd.DataFrame) -> pd.Series:
    """
    Calculate Coefficient of Variation (CV) for each column.

    CV is calculated as (standard deviation / mean) * 100%.
    This function handles potential division by zero by treating the resulting
    NaN values (from 0/0) as 0, as zero standard deviation implies no variance.

    Args:
        data: DataFrame with samples in rows and features in columns.

    Returns:
        A pandas Series containing the CV value for each column.
    """
    mean = data.mean()
    std_dev = data.std()
    cv = (std_dev / mean) * 100
    return cv.fillna(0)

def calculate_rmad(data: pd.DataFrame) -> pd.Series:
    """
    Calculate Relative Median Absolute Deviation (rMAD) for each column.

    rMAD is a robust measure of variability, calculated as:
    (median absolute deviation / median) * 100%.
    It handles potential division by zero similarly to calculate_cv.

    Args:
        data: DataFrame with samples in rows and features in columns.

    Returns:
        A pandas Series containing the rMAD value for each column.
    """
    median_val = data.median()
    mad = (np.abs(data - median_val)).median()
    rmad = (mad / median_val) * 100
    return rmad.fillna(0)

def get_cluster_columns(data: pd.DataFrame) -> List[str]:
    """
    Identify metabolite cluster columns (prefixed with 'N_Cluster' or 'P_Cluster').

    Args:
        data: The input DataFrame.

    Returns:
        A list of column names corresponding to metabolite features.
    """
    return [col for col in data.columns if col.startswith(('N_Cluster', 'P_Cluster'))]

def analyse_variability(data: pd.DataFrame, cluster_columns: List[str]) -> Tuple[float, float]:
    """
    Analyse variability in metabolite data using CV and rMAD metrics.

    Args:
        data: DataFrame containing metabolomics data.
        cluster_columns: A list of metabolite columns to be analysed.

    Returns:
        A tuple containing the percentage of variables exceeding the
        CV and rMAD thresholds, respectively.
    """
    if not cluster_columns:
        return 0.0, 0.0
        
    cluster_data = data[cluster_columns]

    cv_values = calculate_cv(cluster_data)
    rmad_values = calculate_rmad(cluster_data)

    cv_exceeding = (cv_values > CV_THRESHOLD).mean() * 100
    rmad_exceeding = (rmad_values > RMAD_THRESHOLD).mean() * 100

    return cv_exceeding, rmad_exceeding

def process_file(file_path: str) -> Optional[Tuple[str, float, float]]:
    """
    Load a single data file and perform the variability analysis.

    Args:
        file_path: The full path to the input CSV file.

    Returns:
        A tuple containing (file_name, cv_exceeding_pct, rmad_exceeding_pct),
        or None if an error occurs.
    """
    try:
        file_name = os.path.basename(file_path)
        data = pd.read_csv(file_path)
        
        cluster_columns = get_cluster_columns(data)
        if not cluster_columns:
            print(f"Warning: No metabolite cluster columns found in {file_name}. Skipping analysis.")
            return None

        cv_exceeding, rmad_exceeding = analyse_variability(data, cluster_columns)
        
        return file_name, cv_exceeding, rmad_exceeding

    except FileNotFoundError:
        print(f"Error: The file was not found at {file_path}")
        return None
    except pd.errors.EmptyDataError:
        print(f"Error: The file {file_name} is empty. Skipping analysis.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while processing {file_path}: {e}")
        return None

def main():
    """
    Main function to execute the data quality assessment.
    
    It iterates through a list of file paths, processes each file,
    and prints a summary of the variability analysis.
    """
    print("Starting Metabolomics Data Quality Assessment...")
    
    for file_path in FILE_PATHS:
        result = process_file(file_path)
        if result:
            file_name, cv_exceeding, rmad_exceeding = result
            print(f"\n--- Results for {file_name} ---")
            print(f"Percentage of metabolites with CV > {CV_THRESHOLD}%: "
                  f"{cv_exceeding:.2f}%")
            print(f"Percentage of metabolites with rMAD > {RMAD_THRESHOLD}%: "
                  f"{rmad_exceeding:.2f}%")

    print("\nQuality assessment complete.")

if __name__ == "__main__":
    main()
