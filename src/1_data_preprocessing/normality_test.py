"""
Normality Testing for Metabolomics Data
-------------------------------------

This script performs comprehensive normality testing on metabolomics data using two
complementary statistical tests: the Shapiro-Wilk test and the Anderson-Darling test.

It processes multiple CSV files, each representing a different data transformation,
and for each file, it performs the normality tests on all columns containing 'N_Cluster'
in their names.

1.  **Shapiro-Wilk Test**:
    - Ideal for small to medium-sized samples (n < 5000).
    - Tests the null hypothesis that the data was drawn from a normal distribution.

2.  **Anderson-Darling Test**:
    - More sensitive to deviations in the tails of the distribution.
    - Provides a qualitative result (normal or not) based on a 5% significance level.

Input:
- A list of CSV files containing metabolomics data.

Output:
- For each input file, a corresponding CSV file is generated with the normality
  test results. This output includes test statistics and p-values for both methods.
- A warning is printed if the number of variables processed does not match the
  expected count.
"""
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Any
import pandas as pd
from scipy.stats import shapiro, anderson

# --- Configuration ---
CONFIG = {
    'input_dir': Path(r"C:\Users\ms\Desktop\data_chem\data\transformation"),
    'expected_variable_count': 807,
    'significance_level': 5.0,  # For Anderson-Darling test
    'cluster_prefix': 'N_Cluster'
}

def load_data(file_path: Path) -> pd.DataFrame:
    """Loads data from a CSV file."""
    return pd.read_csv(file_path)


def extract_n_cluster_data(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Extracts columns with names containing a given prefix."""
    return df.filter(regex=prefix)


def shapiro_wilk_test(data: pd.DataFrame) -> List[Tuple[str, float, float]]:
    """
    Performs the Shapiro-Wilk normality test on each column of the dataframe.

    Args:
        data: A pandas DataFrame with variables in columns.

    Returns:
        A list of tuples, each containing the column name, Shapiro-Wilk
        statistic, and the corresponding p-value.
    """
    results = []
    for col in data.columns:
        stat, p = shapiro(data[col].dropna())
        results.append((col, stat, p))
    return results


def anderson_darling_test(data: pd.DataFrame, significance_level: float) -> List[Tuple[str, float, bool]]:
    """
    Performs the Anderson-Darling test for normality on each column.

    This implementation provides a qualitative result based on the specified significance level.

    Args:
        data: A pandas DataFrame with variables in columns.
        significance_level: The significance level for the test (e.g., 5.0 for 5%).

    Returns:
        A list of tuples, each containing the column name, Anderson-Darling
        statistic, and a boolean indicating if the data is normal (True) or not (False).
    """
    results = []
    for col in data.columns:
        result = anderson(data[col].dropna())
        statistic = result.statistic
        is_normal = True

        try:
            crit_val_idx = result.significance_level.tolist().index(significance_level)
            critical_value = result.critical_values[crit_val_idx]

            if statistic > critical_value:
                is_normal = False
        except ValueError:
            is_normal = False
        
        results.append((col, statistic, is_normal))
    return results


def save_results_to_csv(results: List, file_path: Path):
    """
    Saves the combined normality test results to a CSV file.

    Args:
        results: A list of lists containing the combined test results.
        file_path: The path for the output CSV file.
    """
    df = pd.DataFrame(results, columns=[
        'Variable',
        'Shapiro_Statistic',
        'Shapiro_p_value',
        'Anderson_Statistic',
        'Anderson_Is_Normal_5_Percent'
    ])
    df.to_csv(file_path, index=False)
    print(f"Normality test results saved to {file_path}")


def process_file(file_path: Path, config: Dict[str, Any]):
    """
    Loads a single data file, performs normality tests, and saves the results.

    Args:
        file_path: Path to the input CSV file.
        config: A dictionary containing configuration parameters.
    """
    print(f"Processing {file_path}...")
    data = load_data(file_path)
    n_cluster_data = extract_n_cluster_data(data, config['cluster_prefix'])

    # Perform normality tests
    shapiro_results = shapiro_wilk_test(n_cluster_data)
    anderson_results = anderson_darling_test(n_cluster_data, config['significance_level'])

    # Combine results
    combined_results = []
    for (var, s_stat, s_p), (_, a_stat, a_is_normal) in zip(shapiro_results, anderson_results):
        combined_results.append([var, s_stat, s_p, a_stat, a_is_normal])

    # Save results
    output_file = file_path.with_name(f"{file_path.stem}_normality_results.csv")
    save_results_to_csv(combined_results, output_file)

    # Validate number of variables
    num_variables = len(n_cluster_data.columns)
    if num_variables != config['expected_variable_count']:
        print(f"Warning: Processed {num_variables} variables, but expected "
              f"{config['expected_variable_count']} in {file_path}")


def main(config: Dict[str, Any]):
    """
    Main function to execute the normality testing pipeline on all specified files.
    """
    file_paths = list(config['input_dir'].glob('*.csv'))
    
    for file_path in file_paths:
        if "normality_results" not in file_path.name:
            if file_path.exists():
                process_file(file_path, config)
            else:
                print(f"Error: File not found at {file_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform normality testing on metabolomics data.")
    parser.add_argument("--input_dir", type=Path, default=CONFIG['input_dir'],
                        help="Directory containing the input CSV files.")
    parser.add_argument("--prefix", type=str, default=CONFIG['cluster_prefix'],
                        help="Prefix for cluster columns.")
    args = parser.parse_args()

    current_config = CONFIG.copy()
    current_config['input_dir'] = args.input_dir
    current_config['cluster_prefix'] = args.prefix
    
    main(current_config)
