"""
Advanced Missing Value Imputation for Metabolomics Data
----------------------------------------------------

This script implements multiple sophisticated imputation methods for handling missing
values in metabolomics data, particularly focusing on P_Cluster measurements.

The implemented methods are:
1.  k-Nearest Neighbors (kNN)
2.  Expectation-Maximization (EM)
3.  Singular Value Decomposition (SVD)
4.  Gaussian Process Regression (GPR)

The script processes a given dataset by first applying one-hot encoding to
categorical variables. Then, it applies each imputation method to the metabolomics
features (P_Cluster columns) and saves the results of each method to a separate
CSV file.

Input:
- A CSV file containing metabolomics data with categorical variables ('Day',
  'Batch', 'Genotype', 'Treatment') and numerical 'P_Cluster' columns.

Output:
- Four CSV files, one for each imputation method, containing the full dataset
  with missing values imputed.
"""
import argparse
from pathlib import Path
import pandas as pd
from sklearn.impute import KNNImputer, IterativeImputer
from fancyimpute import SoftImpute
from sklearn.experimental import enable_iterative_imputer
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from typing import Tuple, List, Dict, Any

# --- Configuration ---
CONFIG = {
    'input_file': Path(r'C:\Users\ms\Desktop\data_chem\data\old_2\p_column_data_r_deducted.csv'),
    'output_base': Path(r'C:\Users\ms\Desktop\data_chem\imputated\p_column_data_r_deducted_imputed_')
}

def load_and_prepare_data(file_path: Path) -> Tuple[pd.DataFrame, List[str]]:
    """
    Load and prepare data for imputation.

    This function loads a CSV file, applies one-hot encoding to categorical
    features, and identifies metabolomics features for imputation.

    Args:
        file_path: The path to the input CSV file.

    Returns:
        A tuple containing:
        - The transformed DataFrame with encoded categorical variables.
        - A list of column names to be imputed.
    """
    data = pd.read_csv(file_path)
    categorical_features = ['Day', 'Batch', 'Genotype', 'Treatment']

    transformer = ColumnTransformer(
        transformers=[
            ('cat', OneHotEncoder(handle_unknown='ignore'), categorical_features)
        ],
        remainder='passthrough'
    )

    data_transformed = transformer.fit_transform(data)
    transformed_cols = transformer.get_feature_names_out()
    transformed_df = pd.DataFrame(data_transformed, columns=transformed_cols)
    
    cols_to_impute = [
        col for col in transformed_cols if 'remainder__P_Cluster' in col
    ]
    
    return transformed_df, cols_to_impute


def impute_and_save(
    data_to_impute: pd.DataFrame,
    columns_to_impute: List[str],
    imputer: Any,
    output_file: Path,
):
    """
    Perform imputation on a copy of the data and save the results.

    Args:
        data_to_impute: A copy of the DataFrame to perform imputation on.
        columns_to_impute: A list of column names to be imputed.
        imputer: An initialized imputer object.
        output_file: The path for saving the imputed data.
    """
    print(f"Starting imputation with {imputer.__class__.__name__}...")
    
    imputed_data = data_to_impute.copy()
    imputed_data[columns_to_impute] = imputer.fit_transform(
        imputed_data[columns_to_impute]
    )

    imputed_data.to_csv(output_file, index=False)
    print(f"Imputed data saved to: {output_file}\n")


def get_imputer_configurations(output_base: Path) -> Dict[str, Dict[str, Any]]:
    """
    Get the configurations for the different imputation methods.

    Args:
        output_base: The base path for the output files.

    Returns:
        A dictionary of imputer configurations.
    """
    return {
        'knn': {
            'imputer': KNNImputer(n_neighbors=5),
            'file': output_base.with_name(f'{output_base.stem}knn.csv')
        },
        'em': {
            'imputer': IterativeImputer(random_state=0),
            'file': output_base.with_name(f'{output_base.stem}em.csv')
        },
        'svd': {
            'imputer': SoftImpute(),
            'file': output_base.with_name(f'{output_base.stem}svd.csv')
        },
        'gpr': {
            'imputer': IterativeImputer(
                estimator=GaussianProcessRegressor(
                    kernel=DotProduct() + WhiteKernel(),
                    alpha=1e-2,
                    random_state=0
                ),
                random_state=0,
                max_iter=15,
                tol=0.001
            ),
            'file': output_base.with_name(f'{output_base.stem}gpr.csv')
        }
    }


def main(config: Dict[str, Path]):
    """
    Main function to run the imputation pipeline.
    """
    prepared_data, cols_to_impute = load_and_prepare_data(config['input_file'])
    imputers = get_imputer_configurations(config['output_base'])

    for method, imputer_config in imputers.items():
        print("-" * 50)
        print(f"Processing imputation method: {method.upper()}")
        
        impute_and_save(
            prepared_data,
            cols_to_impute,
            imputer_config['imputer'],
            imputer_config['file']
        )
    print("-" * 50)
    print("All imputation methods have been applied and saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Advanced missing value imputation for metabolomics data.")
    parser.add_argument("--input_file", type=Path, default=CONFIG['input_file'],
                        help="Path to the input CSV file.")
    parser.add_argument("--output_base", type=Path, default=CONFIG['output_base'],
                        help="Base path for the output files.")
    args = parser.parse_args()

    current_config = {
        'input_file': args.input_file,
        'output_base': args.output_base
    }
    main(current_config)
