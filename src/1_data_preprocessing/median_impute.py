"""
Median-Based Missing Value Imputation for Metabolomics Data
--------------------------------------------------------

This script performs a three-step imputation process for missing values in metabolomics data:
1. Group-wise median imputation based on experimental factors
2. Forward fill for any remaining missing values
3. Backward fill as a final step to ensure complete data

The imputation strategy preserves the data structure by:
- Considering experimental factors (Day, Batch, Genotype, Treatment)
- Using group-specific medians to maintain biological relevance
- Applying sequential imputation methods for comprehensive coverage

Input:
- CSV file containing metabolomics data with P_Cluster columns
- Categorical variables: Day, Batch, Genotype, Treatment
- Missing values represented as NaN

Output:
- CSV file with all missing values imputed
- Original column structure and data types preserved
"""

import sys
from typing import List
import pandas as pd

# --- Configuration ---
# Define grouping columns for imputation
GROUPING_COLUMNS = ['Day', 'Batch', 'Genotype', 'Treatment']
# Define the prefix for metabolite columns
METABOLITE_PREFIX = 'P_Cluster_'


def load_and_prepare_data(file_path: str) -> pd.DataFrame:
    """
    Load the dataset and convert categorical variables to appropriate data types.
    
    Args:
        file_path (str): Path to the input CSV file
        
    Returns:
        pd.DataFrame: Loaded dataframe with properly typed categorical columns
    """
    try:
        data = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{file_path}'")
        sys.exit(1)
    
    # Convert experimental factors to categorical type for efficient grouping
    for col in GROUPING_COLUMNS:
        if col in data.columns:
            data[col] = data[col].astype('category')
    
    return data


def identify_metabolite_columns(data: pd.DataFrame) -> List[str]:
    """
    Identify metabolite columns (P_Cluster) for imputation.
    
    Args:
        data (pd.DataFrame): Input dataframe
        
    Returns:
        list: Column names starting with METABOLITE_PREFIX
    """
    return [col for col in data.columns if col.startswith(METABOLITE_PREFIX)]


def impute_missing_values(data: pd.DataFrame, metabolite_columns: List[str]) -> pd.DataFrame:
    """
    Perform three-step imputation for missing values.
    
    Steps:
    1. Group-wise median imputation
    2. Forward fill remaining NaNs
    3. Backward fill any remaining NaNs
    
    Args:
        data (pd.DataFrame): Input dataframe
        metabolite_columns (list): Columns to impute
        
    Returns:
        pd.DataFrame: Dataframe with imputed values
    """
    imputed_data = data.copy()
    
    # Step 1: Group-wise median imputation
    for column in metabolite_columns:
        imputed_data[column] = (imputed_data.groupby(GROUPING_COLUMNS)[column]
                               .transform(lambda x: x.fillna(x.median())))
    
    # Step 2 & 3: Forward and backward fill remaining NaNs across all metabolite columns
    imputed_data[metabolite_columns] = imputed_data[metabolite_columns].ffill()
    imputed_data[metabolite_columns] = imputed_data[metabolite_columns].bfill()
    
    return imputed_data


def save_imputed_data(data: pd.DataFrame, output_path: str) -> None:
    """
    Save the imputed data to a CSV file.
    
    Args:
        data (pd.DataFrame): Imputed dataframe
        output_path (str): Path for the output CSV file
    """
    try:
        data.to_csv(output_path, index=False)
    except IOError:
        print(f"Error: Could not write to output file at '{output_path}'")
        sys.exit(1)


def main() -> None:
    """
    Main function to execute the imputation pipeline.
    """
    # --- File Paths ---
    # Note: These paths are hardcoded as per requirement.
    input_path = r'C:\Users\ms\Desktop\data_chem\data\old_2\p_column_data_l_deducted.csv'
    output_path = r'C:\Users\ms\Desktop\data_chem\imputated\p_column_data_l_deducted_imputed_median2.csv'
    
    # Load and prepare data
    data = load_and_prepare_data(input_path)
    
    # Identify columns for imputation
    metabolite_columns = identify_metabolite_columns(data)
    
    # Perform imputation
    imputed_data = impute_missing_values(data, metabolite_columns)
    
    # Save results
    save_imputed_data(imputed_data, output_path)
    
    print(f"Data imputation complete. File saved to '{output_path}'")


if __name__ == "__main__":
    main()
