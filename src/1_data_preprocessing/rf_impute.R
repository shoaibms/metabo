# ==============================================================================
# Missing Value Imputation for Metabolomics Data
#
# Author: [Your Name/Organization]
# Date: [Current Date]
#
# Description:
# This script performs missing value imputation on metabolomics data using two
# primary methods:
#
# 1. Random Forest (missForest): An accurate but computationally intensive
#    method that works well with mixed data types and non-linear relationships.
#
# 2. Predictive Mean Matching (MICE): A memory-efficient method that
#    preserves the original data distribution. It is implemented here with
#    column-wise chunking to handle high-dimensional datasets.
#
# The script is structured to be run from top to bottom, with all
# configurations set in the initial section.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
# All user-configurable paths and parameters are defined in this section.
# NOTE: This script uses hardcoded absolute paths as requested. For improved
# portability, consider using relative paths or command-line arguments.

CONFIG <- list(
  # Seed for random number generation to ensure reproducibility
  seed = 2024,
  
  # Number of CPU cores for parallel processing
  # Using detectCores() - 1 is a safe default for most systems.
  cores = max(1, parallel::detectCores() - 1),
  
  # Input and output paths for Random Forest imputation
  rf_imputation = list(
    input = "C:/Users/ms/Desktop/data_chem/data/outlier/p_l_rf_outliers_removed.csv",
    output = "C:/Users/ms/Desktop/data_chem/data/outlier/p_l_if.csv",
    cluster_prefix = "P_Cluster"
  ),
  
  # Input and output for first MICE imputation (100-column chunks)
  mice_imputation_1 = list(
    input = "C:/Users/ms/Desktop/data_chem/data/old_2/p_column_data_r_deducted.csv",
    output = "C:/Users/ms/Desktop/data_chem/imputated/p_column_data_r_deducted_pmm.csv",
    cluster_prefix = "P_Cluster",
    chunk_size = 100
  ),
  
  # Input and output for second MICE imputation (200-column chunks)
  mice_imputation_2 = list(
    input = "C:/Users/ms/Desktop/data_chem/imputated/n_column_data_r_deducted_imputed_pmm3b.csv",
    output = "C:/Users/ms/Desktop/data_chem/imputated/n_column_data_r_deducted_imputed_pmm3c.csv",
    cluster_prefix = "N_Cluster",
    chunk_size = 200
  )
)


# ------------------------------------------------------------------------------
# 2. Package Management
# ------------------------------------------------------------------------------
# This section checks for required packages and installs any that are missing.

install_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    message("Installing required packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, repos = "http://cran.rstudio.com/")
  }
}

required_packages <- c("dplyr", "missForest", "mice", "doParallel")
install_packages(required_packages)

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(missForest)
  library(mice)
  library(doParallel)
})


# ------------------------------------------------------------------------------
# 3. Utility Functions
# ------------------------------------------------------------------------------

#' Prepare Data for Imputation
#'
#' Reads a CSV file, converts specified columns to factors, and selects
#' columns for imputation based on a cluster prefix.
#'
#' @param file_path Path to the input CSV file.
#' @param cluster_prefix Prefix for identifying cluster columns (e.g., "P_Cluster").
#' @return A list containing the full dataset and the subset for imputation.
prepare_data <- function(file_path, cluster_prefix) {
  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path)
  }
  
  data <- read.csv(file_path)
  
  categorical_vars <- c("Day", "Batch", "Genotype", "Treatment")
  data[categorical_vars] <- lapply(data[categorical_vars], as.factor)
  
  # Select columns containing the cluster prefix, plus the categorical variables
  imputation_subset <- dplyr::select(data, dplyr::contains(cluster_prefix), dplyr::all_of(categorical_vars))
  
  return(list(
    full_data = data,
    imputation_data = imputation_subset
  ))
}

#' Save Imputed Data
#'
#' Saves a data frame to a specified CSV file path.
#'
#' @param data The data frame to save.
#' @param file_path The path for the output CSV file.
save_imputed_data <- function(data, file_path) {
  # Ensure the output directory exists
  output_dir <- dirname(file_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  write.csv(data, file_path, row.names = FALSE)
}


# ------------------------------------------------------------------------------
# 4. Imputation Functions
# ------------------------------------------------------------------------------

#' Perform Random Forest Imputation
#'
#' Imputes missing values using the missForest package with parallel processing.
#'
#' @param input_path Path to the input data file.
#' @param output_path Path to save the imputed data.
#' @param cluster_prefix The column prefix for selecting data to impute.
#' @param n_cores The number of CPU cores to use for parallelization.
perform_rf_imputation <- function(input_path, output_path, cluster_prefix, n_cores) {
  message("Starting Random Forest imputation for: ", basename(input_path))
  
  data_list <- prepare_data(input_path, cluster_prefix)
  
  registerDoParallel(cores = n_cores)
  message("Parallel processing registered on ", n_cores, " cores.")
  
  imputed_rf <- missForest::missForest(
    data_list$imputation_data,
    maxiter = 5,
    ntree = 50,
    parallelize = 'variables'
  )
  
  stopImplicitCluster()
  
  # Update the original dataframe with the imputed values
  data_list$full_data[, names(data_list$imputation_data)] <- imputed_rf$ximp
  
  save_imputed_data(data_list$full_data, output_path)
  message("RF imputation complete. Saved to: ", output_path)
}


#' Perform PMM Imputation using MICE with Chunking
#'
#' Imputes missing values by processing the data in column-wise chunks using
#' the MICE package with Predictive Mean Matching (PMM).
#'
#' @param input_path Path to the input data file.
#' @param output_path Path to save the imputed data.
#' @param cluster_prefix The column prefix for selecting data.
#' @param chunk_size The number of columns to process in each chunk.
perform_mice_imputation <- function(input_path, output_path, cluster_prefix, chunk_size) {
  message("Starting MICE (PMM) imputation for: ", basename(input_path))
  message("Using a chunk size of ", chunk_size, " columns.")
  
  data_list <- prepare_data(input_path, cluster_prefix)
  imputation_data <- data_list$imputation_data
  
  n_chunks <- ceiling(ncol(imputation_data) / chunk_size)
  
  for (i in 1:n_chunks) {
    message(sprintf("Processing chunk %d/%d...", i, n_chunks))
    
    start_col <- (i - 1) * chunk_size + 1
    end_col <- min(i * chunk_size, ncol(imputation_data))
    cols_to_impute <- start_col:end_col
    
    chunk_data <- imputation_data[, cols_to_impute, drop = FALSE]
    
    # Skip chunk if there are no missing values
    if (!any(is.na(chunk_data))) {
        message("Skipping chunk ", i, " - no missing values found.")
        next
    }

    imputed_mice <- mice::mice(
      chunk_data,
      method = 'pmm',
      m = 5,
      maxit = 5,
      print = FALSE # Set to TRUE for verbose MICE output
    )
    
    completed_chunk <- mice::complete(imputed_mice, action = 1)
    
    # Update the corresponding columns in the main imputation dataset
    imputation_data[, cols_to_impute] <- completed_chunk
  }
  
  # Update the original full dataframe
  data_list$full_data[, names(imputation_data)] <- imputation_data
  
  save_imputed_data(data_list$full_data, output_path)
  message("MICE imputation complete. Saved to: ", output_path)
}


# ------------------------------------------------------------------------------
# 5. Main Execution
# ------------------------------------------------------------------------------

main <- function() {
  message("Starting imputation script.")
  
  # Set the seed for reproducibility of all random processes
  set.seed(CONFIG$seed)
  message("Global seed set to: ", CONFIG$seed)
  
  # --- Run Random Forest Imputation ---
  rf_params <- CONFIG$rf_imputation
  perform_rf_imputation(
    input_path = rf_params$input,
    output_path = rf_params$output,
    cluster_prefix = rf_params$cluster_prefix,
    n_cores = CONFIG$cores
  )
  
  # --- Run First MICE Imputation ---
  mice_params1 <- CONFIG$mice_imputation_1
  perform_mice_imputation(
    input_path = mice_params1$input,
    output_path = mice_params1$output,
    cluster_prefix = mice_params1$cluster_prefix,
    chunk_size = mice_params1$chunk_size
  )
  
  # --- Run Second MICE Imputation ---
  mice_params2 <- CONFIG$mice_imputation_2
  perform_mice_imputation(
    input_path = mice_params2$input,
    output_path = mice_params2$output,
    cluster_prefix = mice_params2$cluster_prefix,
    chunk_size = mice_params2$chunk_size
  )
  
  message("Script finished successfully.")
}

# This ensures the main function runs only when the script is executed directly
if (sys.nframe() == 0) {
  main()
}
