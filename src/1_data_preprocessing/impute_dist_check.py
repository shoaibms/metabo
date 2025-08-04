"""
Imputation Validation Analysis for Metabolomics Data
-------------------------------------------------

This script evaluates the quality of metabolomics data imputation by comparing
original and imputed distributions using multiple statistical approaches:

1. Earth Mover's Distance (Wasserstein): Measures the minimum "work" needed to
   transform one distribution into another.
2. Hellinger Distance: Quantifies the similarity between probability distributions.
3. Visual comparisons using Q-Q plots, ECDFs, and kernel density estimates.

The script generates both numerical metrics and visualization plots to assess
imputation quality comprehensively.

Output:
- CSV file with distance metrics.
- Composite figure with three plots:
  * Q-Q plot comparing quantiles.
  * Empirical Cumulative Distribution Function (ECDF).
  * Kernel Density Estimation (KDE).
"""

import pandas as pd
import numpy as np
from scipy.stats import wasserstein_distance
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Tuple


def hellinger(p: pd.Series, q: pd.Series) -> float:
    """
    Calculate Hellinger distance between two probability distributions.

    The Hellinger distance is bounded between 0 (identical distributions) and
    1 (completely different distributions).

    Args:
        p (pd.Series): Input probability distribution.
        q (pd.Series): Input probability distribution.

    Returns:
        float: Hellinger distance between p and q.
    """
    return np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) / np.sqrt(2)


def load_and_prepare_data(
    original_path: Path, imputed_path: Path
) -> Tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """
    Load and prepare original and imputed datasets for comparison.

    Args:
        original_path (Path): Path to original data CSV.
        imputed_path (Path): Path to imputed data CSV.

    Returns:
        A tuple containing:
        - pd.Series: Flattened original data (mean across N_Cluster columns).
        - pd.Series: Flattened imputed data (mean across N_Cluster columns).
        - pd.Series: Normalized original data for probability comparison.
        - pd.Series: Normalized imputed data for probability comparison.

    Raises:
        FileNotFoundError: If either of the input files does not exist.
    """
    try:
        original_data = pd.read_csv(original_path)
        imputed_data = pd.read_csv(imputed_path)
    except FileNotFoundError as e:
        print(f"Error loading data: {e}")
        raise

    # This assumes that columns of interest contain 'N_Cluster'.
    # This might need to be adjusted for different datasets.
    original_n_cluster = original_data.filter(like='N_Cluster')
    imputed_n_cluster = imputed_data.filter(like='N_Cluster')

    # Calculate mean across clusters
    original_flat = original_n_cluster.mean(axis=1)
    imputed_flat = imputed_n_cluster.mean(axis=1)

    # Normalize to create a probability distribution for Hellinger distance
    original_norm = original_flat / original_flat.sum()
    imputed_norm = imputed_flat / imputed_flat.sum()

    return original_flat, imputed_flat, original_norm, imputed_norm


def calculate_distances(
    original_flat: pd.Series,
    imputed_flat: pd.Series,
    original_norm: pd.Series,
    imputed_norm: pd.Series,
) -> Tuple[float, float]:
    """
    Calculate distance metrics between original and imputed distributions.

    Args:
        original_flat (pd.Series): Flattened original data.
        imputed_flat (pd.Series): Flattened imputed data.
        original_norm (pd.Series): Normalized original distribution.
        imputed_norm (pd.Series): Normalized imputed distribution.

    Returns:
        A tuple containing the Earth Mover's Distance and Hellinger distance.
    """
    emd = wasserstein_distance(original_flat.dropna(), imputed_flat.dropna())
    hellinger_dist = hellinger(original_norm.dropna(), imputed_norm.dropna())
    return emd, hellinger_dist


def create_comparison_plots(
    original_flat: pd.Series, imputed_flat: pd.Series, cmap_name: str = "BuGn"
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create three comparison plots: Q-Q plot, ECDF, and KDE.

    Args:
        original_flat (pd.Series): Data to compare.
        imputed_flat (pd.Series): Data to compare.
        cmap_name (str): Name of colormap to use.

    Returns:
        A tuple containing the Figure and axes objects.
    """
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    cmap = plt.get_cmap(cmap_name)

    # Q-Q Plot
    quantiles_to_plot = np.linspace(0, 1, num=len(original_flat.dropna()))
    original_quantiles = np.quantile(original_flat.dropna(), quantiles_to_plot)
    imputed_quantiles = np.quantile(imputed_flat.dropna(), quantiles_to_plot)

    axs[0].scatter(original_quantiles, imputed_quantiles, edgecolor='k', c=cmap(0.6))
    axs[0].plot(
        [original_quantiles.min(), original_quantiles.max()],
        [original_quantiles.min(), original_quantiles.max()],
        'r--',
    )
    axs[0].set_xlabel('Quantiles of Original Data (without missing value)')
    axs[0].set_ylabel('Quantiles of Imputed Data')
    axs[0].set_title('Q-Q Plot')

    # ECDF Plot
    sns.ecdfplot(original_flat.dropna(), ax=axs[1], label='Original', color=cmap(0.6))
    sns.ecdfplot(imputed_flat.dropna(), ax=axs[1], label='Imputed', color=cmap(0.8))
    axs[1].set_title('Empirical Cumulative Distribution Function')
    axs[1].legend()

    # KDE Plot
    sns.kdeplot(original_flat.dropna(), ax=axs[2], label='Original', fill=True, color=cmap(0.6))
    sns.kdeplot(imputed_flat.dropna(), ax=axs[2], label='Imputed', fill=True, color=cmap(0.8))
    axs[2].set_title('Kernel Density Estimate')
    axs[2].legend()

    plt.tight_layout()
    return fig, axs


def main() -> None:
    """
    Main function to run the imputation validation analysis.

    It defines the file paths, loads the data, calculates distance metrics,
    saves the metrics to a CSV file, and generates and saves comparison plots.
    """
    # --- Configuration ---
    # File paths are hardcoded as per user request for reproducibility.
    # For a more general-purpose script, consider using command-line arguments.
    original_path = Path(r'C:\Users\ms\Desktop\data_chem\data\old_2\n_column_data_r_deducted.csv')
    imputed_path = Path(r'C:\Users\ms\Desktop\data_chem\imputated\complete imputed files\n_column_data_r_deducted_imputed_rf2.csv')
    output_base_path = Path(r'C:\Users\ms\Desktop\data_chem\imputated\impute_validatiion\n_r_impute_validation_pmm3c')

    # Ensure output directory exists to prevent errors on saving.
    output_base_path.parent.mkdir(parents=True, exist_ok=True)

    # --- Analysis ---
    try:
        # Load and prepare data
        original_flat, imputed_flat, original_norm, imputed_norm = load_and_prepare_data(
            original_path, imputed_path
        )

        # Calculate distance metrics
        emd, hellinger_dist = calculate_distances(
            original_flat, imputed_flat, original_norm, imputed_norm
        )

        # Save distance metrics
        metrics_df = pd.DataFrame({
            "Earth Mover's Distance": [emd],
            'Hellinger Distance': [hellinger_dist]
        })
        csv_output_path = output_base_path.with_suffix('.csv')
        metrics_df.to_csv(csv_output_path, index=False)
        print(f"Distance metrics saved to: {csv_output_path}")

        # Create and save visualization
        fig, _ = create_comparison_plots(original_flat, imputed_flat)
        plot_output_path = output_base_path.with_suffix('.png')
        fig.savefig(plot_output_path)
        print(f"Comparison plot saved to: {plot_output_path}")
        plt.show()

    except FileNotFoundError:
        print("\nExecution stopped: One or more data files were not found.")
        print("Please check the file paths in the 'main' function.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")


if __name__ == "__main__":
    main()
