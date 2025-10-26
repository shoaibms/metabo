"""
Metabolomics Data Filtering Script
------------------------------------

This script filters metabolomics data to remove features with poor replication
in control samples, ensuring higher data quality for downstream analysis.

The filtering logic is as follows:
1.  It considers only control samples (where the 'TMT' column is 0).
2.  For each feature (columns starting with 'N_Cluster'), it counts the number
    of missing values within each 'Entry' group.
3.  A feature column is retained only if all 'Entry' groups have fewer than 3
    missing values.
4.  Essential metadata columns are preserved in the final filtered dataset.
"""

import pandas as pd
from typing import List

# A constant list of metadata columns to preserve in the final dataset.
# These columns are essential for context and are not subject to filtering.
METADATA_COLUMNS = [
    'Row_names', 'Vac_id', 'Entry', 'Maximum_RT_Shift_Neg',
    'Tissue', 'Batch', 'TMT', 'Rep', 'Day'
]

def load_metabolomics_data(file_path: str) -> pd.DataFrame:
    """Loads metabolomics data from an Excel file."""
    return pd.read_excel(file_path, engine='openpyxl')

def identify_metabolite_columns(df: pd.DataFrame) -> List[str]:
    """Identifies columns containing metabolite data based on a prefix."""
    return [col for col in df.columns if col.startswith('N_Cluster')]

def filter_columns_by_replication(df: pd.DataFrame, columns_to_check: List[str]) -> List[str]:
    """
    Filters columns based on replication quality in control samples.

    A column is retained if, for all unique 'Entry' values in the control
    samples (TMT=0), the number of missing values is less than 3. This
    is performed efficiently using vectorized pandas operations.

    Args:
        df: The input DataFrame with metabolomics data.
        columns_to_check: A list of column names to evaluate.

    Returns:
        A list of column names that meet the replication criteria.
    """
    if not columns_to_check:
        return []
    
    control_samples = df[df['TMT'] == 0]

    # Group control samples by 'Entry' and count missing values for each column.
    missing_counts = control_samples.groupby('Entry')[columns_to_check].apply(
        lambda x: x.isnull().sum()
    )

    # A column should be dropped if any entry has 3 or more missing values.
    # This creates a boolean mask to identify such columns.
    columns_to_drop_mask = (missing_counts >= 3).any(axis=0)
    
    # Invert the mask to get the columns to keep and return them as a list.
    columns_to_keep = columns_to_drop_mask[~columns_to_drop_mask].index.tolist()
    
    return columns_to_keep

def save_filtered_data(filtered_df: pd.DataFrame, output_path: str) -> None:
    """Saves the filtered DataFrame to a CSV file."""
    filtered_df.to_csv(output_path, index=False)

def main():
    """Main function to execute the data filtering pipeline."""
    # --- Configuration ---
    # User-specific paths are kept as hardcoded values as requested.
    INPUT_PATH = 'C:\\Users\\ms\\Desktop\\data_chem\\combine\\n_column_data_r.xlsx'
    OUTPUT_PATH = 'C:\\Users\\ms\\Desktop\\data_chem\\combine\\n_column_data_r_clean.csv'

    # --- Execution ---
    print(f"Loading data from {INPUT_PATH}...")
    df = load_metabolomics_data(INPUT_PATH)

    print("Identifying metabolite columns to evaluate...")
    metabolite_columns = identify_metabolite_columns(df)
    print(f"Found {len(metabolite_columns)} metabolite columns.")

    print("Filtering columns based on replication criteria in control samples...")
    columns_to_keep = filter_columns_by_replication(df, metabolite_columns)
    
    final_columns = METADATA_COLUMNS + columns_to_keep
    filtered_df = df[final_columns]

    print(f"Saving filtered data to {OUTPUT_PATH}...")
    save_filtered_data(filtered_df, OUTPUT_PATH)
    
    print("\n--- Filtering Summary ---")
    print(f"Kept {len(columns_to_keep)} of {len(metabolite_columns)} metabolite columns.")
    dropped_count = len(metabolite_columns) - len(columns_to_keep)
    print(f"Dropped {dropped_count} columns.")
    print("-------------------------")
    print("Data filtering completed successfully.")


if __name__ == "__main__":
    main()
