"""Metabolomics Data Transformation Evaluation.

This script evaluates the effectiveness of various data transformations by comparing
statistical metrics and visualizations before and after transformation. It generates
multiple plots and calculates key metrics:

1.  Coefficient of Variation (CV)
2.  MA plots (intensity-dependent ratio plots)
3.  Relative Standard Deviation (RSD)
4.  Relative Median Absolute Deviation (rMAD)

The script saves all generated plots to the 'figure/' directory and outputs a
CSV file with the calculated metrics.
"""

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

# --- Constants ---

# File Paths
# The user requested to keep these paths hardcoded.
# For better portability, consider using relative paths or a configuration file.
BASE_DATA_PATH = Path(r"C:\Users\ms\Desktop\data_chem\data\transformation")
OUTPUT_FIGURE_DIR = Path("figure/transform_metrics")

FILE_PATHS = {
    "original": BASE_DATA_PATH / "n_l_if.csv",
    "anscombe": BASE_DATA_PATH / "n_l_if_anscombe.csv",
    "asinh": BASE_DATA_PATH / "n_l_if_asinh.csv",
    "boxcox": BASE_DATA_PATH / "n_l_if_boxcox.csv",
    "glog": BASE_DATA_PATH / "n_l_if_glog.csv",
    "log": BASE_DATA_PATH / "n_l_if_log.csv",
    "sqrt": BASE_DATA_PATH / "n_l_if_sqrt.csv",
    "yeojohnson": BASE_DATA_PATH / "n_l_if_yeojohnson.csv",
}

# Data extraction
N_CLUSTER_REGEX = "N_Cluster"

# Plotting Styles
CMAP = "BuGn"
TITLE_SIZE = 16
LABEL_SIZE = 14
TICK_SIZE = 12
BEFORE_COLOR = "#8ae391"
AFTER_COLOR = "#0c3b10"


# --- Data Loading ---


def load_data(file_path: Path) -> pd.DataFrame:
    """Load data from a CSV file.

    Args:
        file_path: Path to the CSV file.

    Returns:
        A pandas DataFrame with the loaded data.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Data file not found at: {file_path}")
    return pd.read_csv(file_path)


def extract_n_cluster_data(df: pd.DataFrame) -> pd.DataFrame:
    """Extract N_Cluster columns using a regex pattern."""
    return df.filter(regex=N_CLUSTER_REGEX)


# --- Metric Calculation ---


def calculate_rsd(data: pd.DataFrame) -> np.ndarray:
    """Calculate Relative Standard Deviation (RSD).

    RSD = (standard deviation / mean) * 100
    """
    mean = np.mean(data, axis=0)
    std_dev = np.std(data, axis=0)
    # Avoid division by zero
    rsd = np.divide(std_dev, mean, out=np.zeros_like(std_dev), where=mean != 0) * 100
    return rsd


def calculate_rmad(data: pd.DataFrame) -> np.ndarray:
    """Calculate Relative Median Absolute Deviation (rMAD).

    rMAD = (median(|x - median(x)|) / median(x)) * 100

    This metric is more robust to outliers than RSD.
    """
    median = np.median(data, axis=0)
    mad = np.median(np.abs(data - median), axis=0)
    # Avoid division by zero
    rmad = np.divide(mad, median, out=np.zeros_like(mad), where=median != 0) * 100
    return rmad


# --- Plotting Functions ---


def plot_cv(
    before_df: pd.DataFrame, after_df: pd.DataFrame, title: str, output_dir: Path
) -> Tuple[pd.Series, pd.Series]:
    """Plot Coefficient of Variation comparison.

    CV = (standard deviation / mean) represents relative variability.
    Lower values indicate better precision.
    """
    cv_before = before_df.std() / before_df.mean()
    cv_after = after_df.std() / after_df.mean()

    plt.figure(figsize=(7, 6))
    plt.scatter(cv_before, cv_after, alpha=0.5, c=cv_after, cmap=CMAP)

    max_val = max(cv_before.max(), cv_after.max())
    plt.plot([0, max_val], [0, max_val], "r--")

    plt.title(f"CV Before vs After - {title}", fontsize=TITLE_SIZE)
    plt.xlabel("CV Before Transformation", fontsize=LABEL_SIZE)
    plt.ylabel("CV After Transformation", fontsize=LABEL_SIZE)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    plt.grid(True)
    plt.colorbar(label="CV After")
    plt.savefig(output_dir / f"cv_comparison_{title.lower().replace(' ', '_')}.png")
    plt.close()

    return cv_before, cv_after


def plot_ma_transform(
    before_df: pd.DataFrame, after_df: pd.DataFrame, title: str, output_dir: Path
) -> None:
    """Generate MA plots (intensity-dependent ratio plots).

    MA plots visualize intensity-dependent bias using:
    - A: average of log intensities
    - M: log ratio of intensities
    This implementation uses the first two columns for the plot.
    """
    # Add 1 to avoid log(0)
    before = before_df + 1
    after = after_df + 1

    M_before = np.log2(before.iloc[:, 0]) - np.log2(before.iloc[:, 1])
    A_before = 0.5 * (np.log2(before.iloc[:, 0]) + np.log2(before.iloc[:, 1]))
    M_after = np.log2(after.iloc[:, 0]) - np.log2(after.iloc[:, 1])
    A_after = 0.5 * (np.log2(after.iloc[:, 0]) + np.log2(after.iloc[:, 1]))

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))

    # Before transformation
    scatter1 = axes[0].scatter(A_before, M_before, alpha=0.5, c=M_before, cmap=CMAP)
    lowess_before = sm.nonparametric.lowess(M_before, A_before, frac=0.3)
    axes[0].plot(lowess_before[:, 0], lowess_before[:, 1], color="red")
    axes[0].set_title(f"MA-plot Before - {title}", fontsize=TITLE_SIZE)
    axes[0].set_xlabel("A (average log-intensity)", fontsize=LABEL_SIZE)
    axes[0].set_ylabel("M (log ratio)", fontsize=LABEL_SIZE)
    axes[0].tick_params(axis="both", which="major", labelsize=TICK_SIZE)
    fig.colorbar(scatter1, ax=axes[0], label="M Before")

    # After transformation
    scatter2 = axes[1].scatter(A_after, M_after, alpha=0.5, c=M_after, cmap=CMAP)
    lowess_after = sm.nonparametric.lowess(M_after, A_after, frac=0.3)
    axes[1].plot(lowess_after[:, 0], lowess_after[:, 1], color="red")
    axes[1].set_title(f"MA-plot After - {title}", fontsize=TITLE_SIZE)
    axes[1].set_xlabel("A (average log-intensity)", fontsize=LABEL_SIZE)
    axes[1].set_ylabel("M (log ratio)", fontsize=LABEL_SIZE)
    axes[1].tick_params(axis="both", which="major", labelsize=TICK_SIZE)
    fig.colorbar(scatter2, ax=axes[1], label="M After")

    plt.tight_layout()
    plt.savefig(output_dir / f"ma_plot_{title.lower().replace(' ', '_')}.png")
    plt.close()


def plot_metric_distribution(
    before_df: pd.DataFrame,
    after_df: pd.DataFrame,
    metric_name: str,
    metric_func,
    title: str,
    output_dir: Path,
) -> Tuple[np.ndarray, np.ndarray]:
    """Generate and save boxplot and histogram for a given metric."""
    metric_before = metric_func(before_df)
    metric_after = metric_func(after_df)

    results_df = pd.DataFrame(
        {f"{metric_name}_Before": metric_before, f"{metric_name}_After": metric_after}
    )
    results_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    results_df.dropna(inplace=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sns.boxplot(
        data=[results_df[f"{metric_name}_Before"], results_df[f"{metric_name}_After"]],
        ax=axes[0],
        palette=[BEFORE_COLOR, AFTER_COLOR],
    )
    axes[0].set_xticklabels(["Before", "After"])
    axes[0].set_title(f"Boxplot of {metric_name} - {title}", fontsize=TITLE_SIZE)
    axes[0].set_ylabel(f"{metric_name} (%)", fontsize=LABEL_SIZE)
    axes[0].tick_params(axis="both", which="major", labelsize=TICK_SIZE)

    axes[1].hist(
        [results_df[f"{metric_name}_Before"], results_df[f"{metric_name}_After"]],
        bins=30,
        alpha=0.7,
        label=["Before", "After"],
        color=[BEFORE_COLOR, AFTER_COLOR],
    )
    axes[1].set_title(f"Histogram of {metric_name} - {title}", fontsize=TITLE_SIZE)
    axes[1].set_ylabel("Frequency", fontsize=LABEL_SIZE)
    axes[1].legend(fontsize=LABEL_SIZE)
    axes[1].tick_params(axis="both", which="major", labelsize=TICK_SIZE)

    plt.tight_layout()
    plt.savefig(
        output_dir
        / f"{metric_name.lower()}_dist_{title.lower().replace(' ', '_')}.png"
    )
    plt.close()

    return metric_before, metric_after


# --- Main Execution ---


def save_metrics_to_csv(metrics: Dict[str, List], file_path: Path) -> None:
    """Save computed metrics to a CSV file."""
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(file_path, index=False)
    print(f"Metrics saved to {file_path}")


def main() -> None:
    """Main execution function.

    1. Creates output directory for figures.
    2. Loads original and transformed datasets.
    3. Generates visualization plots for each transformation.
    4. Calculates and collects metrics.
    5. Saves the aggregated metrics to a CSV file.
    """
    OUTPUT_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Saving figures to {OUTPUT_FIGURE_DIR.resolve()}")

    try:
        original_df = load_data(FILE_PATHS["original"])
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    original_n_cluster = extract_n_cluster_data(original_df)

    metrics_summary = []

    for key, path in FILE_PATHS.items():
        if key == "original":
            continue

        print(f"Processing transformation: {key}...")
        try:
            transformed_df = load_data(path)
        except FileNotFoundError as e:
            print(f"Skipping {key}: {e}")
            continue

        transformed_n_cluster = extract_n_cluster_data(transformed_df)
        title = key.replace("_", " ").title()

        cv_before, cv_after = plot_cv(
            original_n_cluster, transformed_n_cluster, title, OUTPUT_FIGURE_DIR
        )
        plot_ma_transform(
            original_n_cluster, transformed_n_cluster, title, OUTPUT_FIGURE_DIR
        )
        rsd_before, rsd_after = plot_metric_distribution(
            original_n_cluster,
            transformed_n_cluster,
            "RSD",
            calculate_rsd,
            title,
            OUTPUT_FIGURE_DIR,
        )
        rmad_before, rmad_after = plot_metric_distribution(
            original_n_cluster,
            transformed_n_cluster,
            "rMAD",
            calculate_rmad,
            title,
            OUTPUT_FIGURE_DIR,
        )

        metrics_summary.append(
            {
                "Transformation": title,
                "Mean_CV_Before": cv_before.mean(),
                "Mean_CV_After": cv_after.mean(),
                "Mean_RSD_Before": np.mean(rsd_before),
                "Mean_RSD_After": np.mean(rsd_after),
                "Mean_rMAD_Before": np.mean(rmad_before),
                "Mean_rMAD_After": np.mean(rmad_after),
            }
        )

    metrics_df = pd.DataFrame(metrics_summary)
    output_csv_path = BASE_DATA_PATH / "transformation_metrics_summary.csv"
    save_metrics_to_csv(metrics_df, output_csv_path)
    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
