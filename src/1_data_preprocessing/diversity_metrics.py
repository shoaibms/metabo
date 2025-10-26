"""Imputation Quality Assessment Using Diversity Metrics
--------------------------------------------------

This script evaluates the quality of different imputation methods (median and Random Forest)
using ecological diversity metrics. These metrics help assess how well the imputation
preserves the data's statistical properties and distribution characteristics.

Metrics calculated:
1. Richness: The number of unique values, which assesses value diversity.
2. Shannon Entropy: Measures the uncertainty and randomness in the distribution.
3. Simpson's Diversity Index: The probability that two randomly selected values are different.
4. Sparsity: The proportion of values equal to the mode, measuring imputation uniformity.

Input:
- Two CSV files containing imputed data from median and Random Forest methods.
- The files should contain 'N_Cluster' columns for analysis.

Output:
- A CSV file with comparative metrics for both imputation methods.
"""

import pandas as pd
import numpy as np
from scipy.stats import entropy
from pathlib import Path
from typing import Tuple, List

CLUSTER_COLUMN_PREFIX = 'N_Cluster'

def load_imputed_datasets(median_path: Path, rf_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and prepare imputed datasets for analysis.

    Args:
        median_path: Path to the median-imputed data.
        rf_path: Path to the RF-imputed data.

    Returns:
        A tuple containing the numeric data from both datasets (median, rf).
    """
    median_data = pd.read_csv(median_path, low_memory=False)
    rf_data = pd.read_csv(rf_path, low_memory=False)
    
    n_cluster_cols: List[str] = [col for col in rf_data.columns if col.startswith(CLUSTER_COLUMN_PREFIX)]
    
    median_numeric = median_data[n_cluster_cols].select_dtypes(include=[np.number])
    rf_numeric = rf_data[n_cluster_cols].select_dtypes(include=[np.number])
    
    return median_numeric, rf_numeric

def calculate_richness(data: pd.DataFrame) -> float:
    """Calculate the average number of unique values across columns.

    Higher richness indicates a more diverse value distribution.
    """
    return data.nunique().mean()

def calculate_shannon_entropy(data: pd.DataFrame) -> float:
    """Calculate Shannon entropy to measure distribution uncertainty.

    Higher entropy suggests a more even distribution of values and better
    preservation of data variability.
    """
    def entropy_of_series(series: pd.Series) -> float:
        counts = series.value_counts()
        return entropy(counts)
    
    return data.apply(entropy_of_series).mean()

def calculate_simpsons_index(data: pd.DataFrame) -> float:
    """Calculate Simpson's diversity index (1 - D).

    This index ranges from 0 (no diversity) to 1 (infinite diversity) and
    represents the probability that two randomly chosen values are different.
    """
    def simpsons_index_of_series(series: pd.Series) -> float:
        counts = series.value_counts()
        n = sum(counts)
        if n == 0:
            return 0.0
        return 1 - sum((count / n) ** 2 for count in counts)

    return data.apply(simpsons_index_of_series).mean()

def calculate_sparsity(imputed_data: pd.DataFrame) -> float:
    """Calculate the average proportion of mode values in each column.

    Lower sparsity suggests better preservation of data variability. High sparsity
    might indicate an over-reliance on central tendencies during imputation.
    """
    modes = imputed_data.mode().iloc[0]
    sparsity = imputed_data.eq(modes, axis=1).mean()
    return sparsity.mean()

def evaluate_imputation_quality(median_data: pd.DataFrame, rf_data: pd.DataFrame) -> pd.DataFrame:
    """Calculate all diversity metrics for both imputation methods.

    Returns:
        A DataFrame comparing all calculated metrics for both methods.
    """
    metrics = {
        'Metric': ['Richness', 'Shannon Entropy', "Simpson's Diversity Index", 'Sparsity'],
        'Median Imputed': [
            calculate_richness(median_data),
            calculate_shannon_entropy(median_data),
            calculate_simpsons_index(median_data),
            calculate_sparsity(median_data)
        ],
        'RF Imputed': [
            calculate_richness(rf_data),
            calculate_shannon_entropy(rf_data),
            calculate_simpsons_index(rf_data),
            calculate_sparsity(rf_data)
        ]
    }
    
    return pd.DataFrame(metrics)

def main() -> None:
    """Main function to execute the imputation quality assessment."""
    # Define input file paths
    median_path = Path(r'C:\Users\ms\Desktop\data_chem\imputated\complete imputed files\n_column_data_r_deducted_imputed_median2.csv')
    rf_path = Path(r'C:\Users\ms\Desktop\data_chem\imputated\complete imputed files\n_column_data_r_deducted_imputed_rf2.csv')
    
    # Load and prepare data
    median_data, rf_data = load_imputed_datasets(median_path, rf_path)
    
    # Calculate and compile all metrics
    results = evaluate_imputation_quality(median_data, rf_data)
    
    # Save results
    output_path = Path(r'C:\Users\ms\Desktop\data_chem\imputated\complete imputed files\imputation_comparison_results.csv')
    results.to_csv(output_path, index=False)
    
    # Print results for immediate review
    print("\nImputation Quality Metrics:")
    print(results)
    print(f"\nResults saved to: {output_path}")

if __name__ == "__main__":
    main()
