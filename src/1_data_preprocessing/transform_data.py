"""
Data Transformation Script for Metabolomics Analysis
-------------------------------------------------

This script applies various data transformations commonly used in metabolomics 
data analysis to address non-normality, heteroscedasticity, and improve 
data distribution properties.

The script performs the following transformations:
1.  Log transformation: ln(x + 1) to handle zeros and reduce right skewness.
2.  Square root: √x to reduce right skewness, less aggressive than log.
3.  Box-Cox: Parametric power transformation family to normalize data.
4.  Yeo-Johnson: Similar to Box-Cox but handles negative values.
5.  Asinh (inverse hyperbolic sine): Similar to log but handles zeros/negative values.
6.  Generalized log (glog): Handles zeros/small values better than regular log.
7.  Anscombe: Variance-stabilizing transformation, especially for count data.

Input:
- A CSV file containing metabolomics data. Columns to be transformed should be 
  prefixed with 'N_Cluster'.
- Non-cluster columns are preserved without change.

Output:
- Separate CSV files for each transformation method, saved in the same directory
  as the input file.
- The original column names and structure are preserved in the output files.
"""

import pandas as pd
import numpy as np
from scipy.stats import boxcox
from sklearn.preprocessing import PowerTransformer
from pathlib import Path
from typing import List, Dict, Callable, Tuple

FILE_PATH = Path(r'C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if.csv')

def load_data(file_path: Path) -> Tuple[pd.DataFrame, List[str]]:
    """
    Loads the dataset from a CSV file and identifies cluster columns.

    Args:
        file_path: The path to the input CSV file.

    Returns:
        A tuple containing the loaded DataFrame and a list of cluster column names.
        
    Raises:
        FileNotFoundError: If the file at file_path does not exist.
    """
    try:
        data = pd.read_csv(file_path)
        cluster_columns = [col for col in data.columns if col.startswith('N_Cluster')]
        print(f"Loaded data with {len(cluster_columns)} cluster columns to transform from {file_path}.")
        return data, cluster_columns
    except FileNotFoundError:
        print(f"Error: Input file not found at {file_path}")
        raise

def apply_log_transform(data: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Applies natural logarithm transformation: ln(x + 1).

    The addition of 1 prevents issues with zero values while preserving order.
    """
    transformed = data.copy()
    transformed[columns] = np.log1p(transformed[columns])
    return transformed

def apply_sqrt_transform(data: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Applies square root transformation: √x.

    This is less aggressive than log transformation for reducing right skewness.
    """
    transformed = data.copy()
    transformed[columns] = np.sqrt(transformed[columns])
    return transformed

def apply_boxcox_transform(data: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Applies Box-Cox transformation to approximate normality.

    A small constant (1) is added to handle zero values before transformation.
    """
    transformed = data.copy()
    for col in columns:
        # Box-Cox requires positive data, adding 1 to handle zeros.
        transformed[col], _ = boxcox(transformed[col] + 1)
    return transformed

def apply_yeojohnson_transform(data: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Applies Yeo-Johnson transformation.

    This is an extension of Box-Cox that can handle zero and negative values.
    """
    transformed = data.copy()
    pt = PowerTransformer(method='yeo-johnson')
    transformed[columns] = pt.fit_transform(transformed[columns])
    return transformed

def apply_asinh_transform(data: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Applies inverse hyperbolic sine (asinh) transformation.

    Similar to log but handles zeros and negative values naturally.
    """
    transformed = data.copy()
    transformed[columns] = np.arcsinh(transformed[columns])
    return transformed

def apply_glog_transform(data: pd.DataFrame, columns: List[str], a: float = 0.75) -> pd.DataFrame:
    """
    Applies generalized logarithm transformation with parameter 'a'.

    This transformation handles small and zero values better than regular log.
    The formula is: log(x + sqrt(x^2 + a^2)).

    Args:
        data: The input DataFrame.
        columns: A list of column names to transform.
        a: A small constant to stabilize variance (default: 0.75).
    """
    transformed = data.copy()
    x = transformed[columns]
    transformed[columns] = np.log(x + np.sqrt(x**2 + a**2))
    return transformed

def apply_anscombe_transform(data: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Applies Anscombe transformation: 2 * √(x + 3/8).

    This is a variance-stabilizing transformation, particularly for count data.
    """
    transformed = data.copy()
    transformed[columns] = 2 * np.sqrt(transformed[columns] + 3/8)
    return transformed

def save_transformed_data(transformed_data: pd.DataFrame, original_path: Path, suffix: str) -> Path:
    """
    Saves transformed data to a new CSV file with a given suffix.
    """
    output_path = original_path.with_name(f"{original_path.stem}_{suffix}.csv")
    transformed_data.to_csv(output_path, index=False)
    print(f"Saved {suffix} transformed data to: {output_path}")
    return output_path

def main() -> None:
    """
    Main function to execute the data transformation pipeline.
    """
    try:
        data, cluster_columns = load_data(FILE_PATH)
    except FileNotFoundError:
        return # Exit if the file is not found
    
    if not cluster_columns:
        print("No columns starting with 'N_Cluster' found. Exiting.")
        return

    # Define the transformations to apply
    transformations: Dict[str, Callable[[pd.DataFrame, List[str]], pd.DataFrame]] = {
        'log': apply_log_transform,
        'sqrt': apply_sqrt_transform,
        'boxcox': apply_boxcox_transform,
        'yeojohnson': apply_yeojohnson_transform,
        'asinh': apply_asinh_transform,
        'glog': apply_glog_transform,
        'anscombe': apply_anscombe_transform
    }
    
    # Apply and save each transformation
    for name, transform_func in transformations.items():
        print(f"Applying {name} transformation...")
        transformed = transform_func(data.copy(), cluster_columns)
        save_transformed_data(transformed, FILE_PATH, name)
    
    print("\nAll transformations completed successfully.")

if __name__ == "__main__":
    main()
