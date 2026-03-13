"""
Little's MCAR (Missing Completely At Random) Test Implementation.

This script implements Little's MCAR test to evaluate whether missing data in
a dataset can be considered Missing Completely At Random (MCAR). This is a
key assumption for many data imputation techniques.

The test's hypotheses are:
- Null hypothesis (H0): The data is Missing Completely At Random.
- Alternative hypothesis (H1): The data is not Missing Completely At Random.

A low p-value (e.g., < 0.05) suggests that the null hypothesis should be
rejected, indicating that the data is not MCAR.

Input:
- A CSV file containing the data to be analyzed. The script specifically
  looks for numeric columns prefixed with 'N_Cluster_'.

Output:
- A CSV file with the p-value from Little's MCAR test.
- A log message to the console with the test result and its interpretation.
"""

import logging
import pandas as pd
import numpy as np
from pyampute.exploration.mcar_statistical_tests import MCARTest
from typing import List

# --- Configuration ---
# Update these paths to match your file system structure.
DATA_PATH = r'C:\Users\USER\Desktop\data_chem\data\n_column_data_l_deducted.csv'
RESULT_PATH = r'C:\Users\USER\Desktop\data_chem\data\littles_mcar_test_results.csv'

# --- Logging Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def load_metabolomics_data(file_path: str) -> pd.DataFrame:
    """
    Loads metabolomics data and extracts numeric columns prefixed with 'N_Cluster_'.

    Args:
        file_path: The path to the CSV file.

    Returns:
        A DataFrame containing only the filtered numeric columns.

    Raises:
        FileNotFoundError: If the specified file_path does not exist.
        ValueError: If no 'N_Cluster_' columns are found.
    """
    logging.info(f"Loading data from {file_path}...")
    try:
        data = pd.read_csv(file_path)
    except FileNotFoundError:
        logging.error(f"Error: The file was not found at {file_path}")
        raise

    n_cluster_columns: List[str] = [
        col for col in data.columns
        if col.startswith('N_Cluster_') and np.issubdtype(data[col].dtype, np.number)
    ]

    if not n_cluster_columns:
        raise ValueError("No numeric 'N_Cluster_' columns found in the dataset.")

    logging.info(f"Found {len(n_cluster_columns)} 'N_Cluster_' columns.")
    return data[n_cluster_columns]


def perform_littles_mcar_test(data: pd.DataFrame) -> float:
    """
    Performs Little's MCAR test on the given DataFrame.

    Args:
        data: The DataFrame to be tested.

    Returns:
        The p-value from Little's MCAR test.
    """
    logging.info("Performing Little's MCAR test...")
    mcar_tester = MCARTest(method="little")
    p_value = mcar_tester.little_mcar_test(data)
    logging.info(f"Little's MCAR test completed. P-value: {p_value:.4f}")
    return p_value


def save_test_results(p_value: float, output_path: str) -> None:
    """
    Saves the MCAR test p-value to a CSV file.

    Args:
        p_value: The p-value from the test.
        output_path: The path to save the results CSV file.
    """
    logging.info(f"Saving results to {output_path}...")
    result_df = pd.DataFrame({'P-Value': [p_value]})
    try:
        result_df.to_csv(output_path, index=False)
        logging.info("Results saved successfully.")
    except IOError as e:
        logging.error(f"Failed to save results to {output_path}: {e}")
        raise


def main() -> None:
    """
    Main function to run the MCAR test pipeline.
    """
    try:
        # Load and preprocess data
        n_cluster_data = load_metabolomics_data(DATA_PATH)

        # Perform Little's MCAR test
        p_value = perform_littles_mcar_test(n_cluster_data)

        # Save results
        save_test_results(p_value, RESULT_PATH)

        # Print interpretation
        logging.info("\n--- Little's MCAR Test Interpretation ---")
        logging.info(f"P-value: {p_value:.4f}")
        if p_value < 0.05:
            logging.warning("The data is likely NOT Missing Completely At Random (MCAR).")
        else:
            logging.info("There is no evidence to suggest the data is not MCAR.")
        logging.info("------------------------------------")

    except FileNotFoundError:
        logging.error("Analysis aborted because a file was not found.")
    except ValueError as e:
        logging.error(f"Analysis aborted due to a value error: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()
