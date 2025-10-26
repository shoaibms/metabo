"""
Comprehensive Outlier Analysis for Metabolomics Data
--------------------------------------------------

This script implements multiple outlier detection methods and their visualization
for metabolomics data analysis. It combines traditional statistical approaches
with modern machine learning techniques for robust outlier detection.

Methods implemented:
1. Z-Score: A traditional statistical method that uses standard deviations.
2. IQR: The interquartile range method for detecting statistical outliers.
3. Isolation Forest: An unsupervised learning method for anomaly detection.
4. Elliptic Envelope: Assumes a Gaussian distribution for outlier detection.
5. Mahalanobis Distance: Accounts for the covariance structure in multivariate data.
6. Robust PCA: Principal Component Analysis with robust covariance estimation.
7. Local Outlier Factor: A density-based outlier detection method.

Visualizations:
- PCA projection with outlier highlighting.
- t-SNE projection with outlier highlighting.
- Performance evaluation metrics for each method.
"""

import logging
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import IsolationForest
from sklearn.covariance import EllipticEnvelope
from sklearn.neighbors import LocalOutlierFactor
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import List, Tuple, Dict
from pathlib import Path

class Config:
    """Configuration class for the outlier detection script."""
    FILE_PATH = Path(r'C:\Users\ms\Desktop\data_chem\data\impute_missing_value\p_l_rf.csv')
    CATEGORICAL_VARS = ['Day', 'Batch', 'Genotype', 'Treatment', 'Replication']
    CONTAMINATION_RATE = 0.05
    Z_SCORE_THRESHOLD = 5.5
    IQR_MULTIPLIER = 7.0
    RANDOM_STATE = 42
    OUTPUT_CSV_PATH = Path('outlier_detection_evaluation_results.csv')
    PCA_N_COMPONENTS = 2
    TSNE_N_COMPONENTS = 2
    LOF_N_NEIGHBORS = 20
    FIGURE_SIZE = (8, 6)
    FONT_SIZE_TITLE = 16
    FONT_SIZE_LABEL = 14
    FONT_SIZE_TICKS = 12
    FONT_SIZE_LEGEND = 12
    
    # Custom color map for visualizations
    NORMAL_COLOR = '#bfdbc1'
    OUTLIER_COLOR = '#425944'

def setup_logging():
    """Configure logging for the script."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_and_prepare_data(
    file_path: Path, categorical_vars: List[str]
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Load and preprocess the metabolomics data.

    Args:
        file_path: The path to the input CSV file.
        categorical_vars: A list of categorical variables to encode.

    Returns:
        A tuple containing:
        - The original DataFrame with encoded categorical variables.
        - Standardized data for 'P_Cluster' columns.
        - A list of 'P_Cluster' column names.
        
    Raises:
        FileNotFoundError: If the file at file_path does not exist.
        KeyError: If a variable in categorical_vars is not in the DataFrame.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"The file '{file_path}' was not found.")
    
    df = pd.read_csv(file_path)

    label_encoders = {}
    for var in categorical_vars:
        if var not in df.columns:
            raise KeyError(f"Categorical variable '{var}' not found in the DataFrame columns.")
        le = LabelEncoder()
        df[var] = le.fit_transform(df[var])
        label_encoders[var] = le
    
    p_cluster_columns = [col for col in df.columns if 'P_Cluster' in col]
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df[p_cluster_columns])
    
    return df, scaled_data, p_cluster_columns

def mahalanobis_distance_outliers(data: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Calculate Mahalanobis distance and identify outliers using an adaptive threshold.

    Args:
        data: The input data matrix (samples x features).

    Returns:
        A tuple containing:
        - A boolean array indicating the outlier status for each sample.
        - The calculated Mahalanobis distances.
        - The outlier detection threshold.
    """
    mean = np.mean(data, axis=0)
    centered_data = data - mean
    
    try:
        cov_matrix = np.cov(centered_data, rowvar=False)
        inv_cov_matrix = np.linalg.inv(cov_matrix)
    except np.linalg.LinAlgError:
        logging.warning("Singular covariance matrix encountered. Using pseudo-inverse.")
        cov_matrix = np.cov(centered_data, rowvar=False)
        inv_cov_matrix = np.linalg.pinv(cov_matrix)

    distances = np.sqrt(np.sum((centered_data @ inv_cov_matrix) * centered_data, axis=1))
    
    median_dist = np.median(distances)
    mad = np.median(np.abs(distances - median_dist))
    threshold = median_dist + 3 * 1.4826 * mad  # 1.4826 scales MAD to be an unbiased estimator of the standard deviation
    outliers = distances > threshold
    
    return outliers, distances, threshold

def robust_pca_outliers(data: np.ndarray, contamination: float, random_state: int, n_components: int) -> np.ndarray:
    """
    Detect outliers using robust PCA with Elliptic Envelope.

    Args:
        data: The input data matrix.
        contamination: The expected proportion of outliers.
        random_state: A seed for reproducibility.
        n_components: The number of principal components to use.

    Returns:
        A boolean array indicating the outlier status for each sample.
    """
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(data)
    robust_cov = EllipticEnvelope(contamination=contamination, random_state=random_state)
    return robust_cov.fit_predict(scores) == -1

def apply_outlier_methods(data: np.ndarray, config: Config) -> Dict[str, np.ndarray]:
    """
    Apply multiple outlier detection methods to the data.

    Args:
        data: The input data matrix.
        config: The configuration object.

    Returns:
        A dictionary containing boolean outlier flags for each method.
    """
    z_scores = np.abs((data - data.mean(axis=0)) / data.std(axis=0))
    z_score_outliers = (z_scores > config.Z_SCORE_THRESHOLD).any(axis=1)

    q1 = np.percentile(data, 25, axis=0)
    q3 = np.percentile(data, 75, axis=0)
    iqr = q3 - q1
    iqr_outliers = ((data < (q1 - config.IQR_MULTIPLIER * iqr)) | (data > (q3 + config.IQR_MULTIPLIER * iqr))).any(axis=1)

    isolation_forest = IsolationForest(contamination=config.CONTAMINATION_RATE, random_state=config.RANDOM_STATE)
    elliptic_envelope = EllipticEnvelope(contamination=config.CONTAMINATION_RATE, random_state=config.RANDOM_STATE)
    lof = LocalOutlierFactor(n_neighbors=config.LOF_N_NEIGHBORS, contamination=config.CONTAMINATION_RATE)
    
    mahalanobis_outliers, _, _ = mahalanobis_distance_outliers(data)

    return {
        'Z-Score': z_score_outliers,
        'IQR': iqr_outliers,
        'Isolation Forest': isolation_forest.fit_predict(data) == -1,
        'Elliptic Envelope': elliptic_envelope.fit_predict(data) == -1,
        'Mahalanobis Distance': mahalanobis_outliers,
        'Robust PCA': robust_pca_outliers(data, config.CONTAMINATION_RATE, config.RANDOM_STATE, config.PCA_N_COMPONENTS),
        'Local Outlier Factor': lof.fit_predict(data) == -1
    }

def evaluate_methods(outlier_results: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    Calculate the percentage of samples identified as outliers by each method.
    
    Args:
        outlier_results: A dictionary of outlier detection results.
        
    Returns:
        A dictionary with the percentage of outliers for each method.
    """
    return {method: np.mean(results) * 100 for method, results in outlier_results.items()}

def plot_dimensionality_reduction(
    data: np.ndarray, outliers: np.ndarray, method_name: str, reduction_type: str, config: Config
):
    """
    Visualize outliers using dimensionality reduction (PCA or t-SNE).

    Args:
        data: The input data matrix.
        outliers: A boolean array indicating outlier status.
        method_name: The name of the outlier detection method.
        reduction_type: The type of dimensionality reduction ('PCA' or 't-SNE').
        config: The configuration object.
    """
    if reduction_type == 'PCA':
        reducer = PCA(n_components=config.PCA_N_COMPONENTS, random_state=config.RANDOM_STATE)
        xlabel, ylabel = 'Principal Component 1', 'Principal Component 2'
    elif reduction_type == 't-SNE':
        reducer = TSNE(n_components=config.TSNE_N_COMPONENTS, random_state=config.RANDOM_STATE)
        xlabel, ylabel = 't-SNE 1', 't-SNE 2'
    else:
        raise ValueError("reduction_type must be either 'PCA' or 't-SNE'")

    reduced_data = reducer.fit_transform(data)
    
    plt.figure(figsize=config.FIGURE_SIZE)
    colors = [config.NORMAL_COLOR, config.OUTLIER_COLOR]
    cmap = mcolors.LinearSegmentedColormap.from_list("Custom", colors, N=2)
    
    plt.scatter(reduced_data[:, 0], reduced_data[:, 1], c=outliers, cmap=cmap, alpha=0.8, s=50)
    
    plt.title(f"{reduction_type} plot of outliers detected by {method_name}", fontsize=config.FONT_SIZE_TITLE)
    plt.xlabel(xlabel, fontsize=config.FONT_SIZE_LABEL)
    plt.ylabel(ylabel, fontsize=config.FONT_SIZE_TICKS)
    plt.tick_params(axis='both', which='major', labelsize=config.FONT_SIZE_TICKS)
    
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Normal Data', markerfacecolor=config.NORMAL_COLOR, markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Outlier', markerfacecolor=config.OUTLIER_COLOR, markersize=10)
    ]
    plt.legend(handles=legend_elements, fontsize=config.FONT_SIZE_LEGEND)
    plt.tight_layout()
    plt.show()

def run_analysis(config: Config):
    """
    Run the full outlier detection and visualization pipeline.

    Args:
        config: The configuration object.
    """
    try:
        df, scaled_data, p_cluster_columns = load_and_prepare_data(config.FILE_PATH, config.CATEGORICAL_VARS)
        data_for_detection = pd.DataFrame(scaled_data, columns=p_cluster_columns)
    except (FileNotFoundError, KeyError) as e:
        logging.error(f"Failed to load or preprocess data: {e}")
        return

    logging.info(f"Data loaded successfully with {data_for_detection.shape[0]} samples and {data_for_detection.shape[1]} features.")

    outlier_results = apply_outlier_methods(data_for_detection.values, config)
    
    _, mahalanobis_dist, mahalanobis_thresh = mahalanobis_distance_outliers(data_for_detection.values)
    logging.info("\n--- Mahalanobis Distance Details ---\n"
                 f"Min Distance: {np.min(mahalanobis_dist):.2f}, Max Distance: {np.max(mahalanobis_dist):.2f}\n"
                 f"Adaptive Threshold: {mahalanobis_thresh:.2f}\n"
                 f"Number of outliers detected: {np.sum(outlier_results['Mahalanobis Distance'])}\n"
                 "------------------------------------\n")

    evaluation_results = evaluate_methods(outlier_results)
    evaluation_df = pd.DataFrame.from_dict(evaluation_results, orient='index', columns=['Percentage_Outliers'])
    evaluation_df = evaluation_df.sort_values(by='Percentage_Outliers', ascending=False)
    
    try:
        evaluation_df.to_csv(config.OUTPUT_CSV_PATH)
        logging.info(f"Evaluation results saved to '{config.OUTPUT_CSV_PATH}'")
    except IOError as e:
        logging.error(f"Error saving evaluation results to CSV: {e}")

    logging.info("\n--- Outlier Detection Evaluation ---\n"
                 f"{evaluation_df}\n"
                 "------------------------------------\n")
    
    logging.info("Generating visualizations for each method...")
    for method, outliers in outlier_results.items():
        if np.sum(outliers) > 0:
            plot_dimensionality_reduction(data_for_detection.values, outliers, method, 'PCA', config)
            plot_dimensionality_reduction(data_for_detection.values, outliers, method, 't-SNE', config)
        else:
            logging.info(f"Skipping visualization for {method} as no outliers were detected.")

    logging.info("Outlier detection process completed.")

def main():
    """Main function to set up and run the outlier detection pipeline."""
    setup_logging()
    config = Config()
    logging.info("Starting the outlier detection process...")
    run_analysis(config)

if __name__ == "__main__":
    main()
