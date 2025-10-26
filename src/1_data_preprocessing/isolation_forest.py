"""
Outlier Detection and Removal Using Isolation Forest
-------------------------------------------------

This script implements outlier detection and removal for metabolomics data using
the Isolation Forest algorithm. Isolation Forest is particularly effective for
high-dimensional data and does not make assumptions about the data distribution.

Algorithm details:
- Uses Isolation Forest to identify outliers in metabolomic features.
- Performs column-wise outlier detection after standardisation.
- Replaces identified outliers with NaN for subsequent imputation.
- Uses the decision function scores to determine the outlier threshold.

Parameters:
- contamination: The expected proportion of outliers (set to 0.05 or 5%).
- threshold: Uses the 5th percentile of decision scores as the cutoff.

Input:
- A CSV file containing metabolomics data with 'P_Cluster' columns.
- Non-cluster columns are preserved unchanged.

Output:
- A CSV file with outliers replaced by NaN values.
- The original column structure and non-outlier values are preserved.
"""
import sys
from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler


def load_data(file_path: str) -> Tuple[pd.DataFrame, List[str]]:
    """
    Load the dataset and identify P_Cluster columns.

    Args:
        file_path (str): Path to the input CSV file.

    Returns:
        tuple: A DataFrame and a list of P_Cluster column names.

    Raises:
        SystemExit: If the file is not found.
    """
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {file_path}")
        sys.exit(1)

    cluster_columns = [col for col in df.columns if 'P_Cluster' in col]
    return df, cluster_columns


def initialise_isolation_forest(
    contamination: float = 0.05, random_state: int = 42
) -> IsolationForest:
    """
    Initialise the Isolation Forest model with specified parameters.

    Args:
        contamination (float): Expected proportion of outliers. Defaults to 0.05.
        random_state (int): Random seed for reproducibility. Defaults to 42.

    Returns:
        IsolationForest: An initialised Isolation Forest model.
    """
    return IsolationForest(contamination=contamination, random_state=random_state)


def process_outliers(
    df: pd.DataFrame,
    columns: List[str],
    model: IsolationForest,
    threshold_percentile: int = 5,
) -> pd.DataFrame:
    """
    Detect and replace outliers in specified columns with NaN.

    This function iterates through each specified column, scales the data,
    fits the Isolation Forest model, and identifies outliers based on the
    decision function scores. Outliers are determined by a percentile-based
    threshold on these scores.

    Args:
        df (pd.DataFrame): The input DataFrame.
        columns (list): A list of column names to process for outliers.
        model (IsolationForest): An initialised Isolation Forest model.
        threshold_percentile (int): The percentile for the decision score
                                     threshold. Defaults to 5.

    Returns:
        pd.DataFrame: A DataFrame with outliers replaced by NaN.
    """
    df_cleaned = df.copy()
    scaler = StandardScaler()

    for col in columns:
        # Reshape and scale the data for a single feature
        data_to_check = df_cleaned[col].dropna()
        if data_to_check.empty:
            continue
            
        scaled_data = scaler.fit_transform(data_to_check.values.reshape(-1, 1))

        # Fit the model and get decision scores for anomaly detection
        model.fit(scaled_data)
        decision_scores = model.decision_function(scaled_data)

        # Determine the outlier threshold from the decision scores
        threshold = np.percentile(decision_scores, threshold_percentile)

        # Identify outliers as values below the threshold
        outlier_indices = data_to_check.index[decision_scores < threshold]

        # Replace outliers with NaN
        df_cleaned.loc[outlier_indices, col] = np.nan

    return df_cleaned


def save_cleaned_data(df: pd.DataFrame, output_path: str) -> None:
    """
    Save the cleaned DataFrame to a CSV file.

    Args:
        df (pd.DataFrame): The cleaned DataFrame to save.
        output_path (str): The path for the output CSV file.

    Raises:
        SystemExit: If the file cannot be saved.
    """
    try:
        df.to_csv(output_path, index=False)
        print(f"Cleaned data saved to: {output_path}")
    except IOError as e:
        print(f"Error: Could not save file to {output_path}. Reason: {e}")
        sys.exit(1)


def main() -> None:
    """Main function to run the outlier detection script."""
    # File paths
    input_path = 'C:/Users/ms/Desktop/data_chem/data/outlier/p_l_rf.csv'
    output_path = 'C:/Users/ms/Desktop/data_chem/data/outlier/p_l_rf_outliers_removed.csv'

    # Load data
    df, cluster_columns = load_data(input_path)

    # Initialise model
    model = initialise_isolation_forest()

    # Process outliers
    df_cleaned = process_outliers(df, cluster_columns, model)

    # Save results
    save_cleaned_data(df_cleaned, output_path)


if __name__ == "__main__":
    main()
