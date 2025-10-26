# Section 1: Fixed Data Processing & Validation
library(tidyverse)
library(stats)

# File paths
base_input_path <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/group/merge2"
base_output_path <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"

# Create output directory if doesn't exist
dir.create(base_output_path, recursive = TRUE, showWarnings = FALSE)

# Modified validation function
validate_metabolite_data <- function(data) {
  validation_results <- list(
    valid = TRUE,
    messages = character()
  )
  
  # Check if data is empty
  if(nrow(data) == 0 || ncol(data) == 0) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages,
                                     "Error: Empty dataset")
    return(validation_results)
  }
  
  # Check for required columns
  required_cols <- c("Metabolite", "Tissue.type")
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if(length(missing_cols) > 0) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages,
                                     paste("Error: Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  return(validation_results)
}

# Modified load and process function
load_process_data <- function() {
  tryCatch({
    # Load data with explicit encoding
    filtered_metabolites <- read.csv(
      file.path(base_input_path, "Merged_VIP_hub_name_only.csv"),
      stringsAsFactors = FALSE,
      check.names = FALSE,
      encoding = "UTF-8"
    )
    
    merged_data <- read.csv(
      file.path(base_input_path, "Merged_VIP_hub.csv"),
      stringsAsFactors = FALSE,
      check.names = FALSE,
      encoding = "UTF-8"
    )
    
    # Validate input data
    filtered_validation <- validate_metabolite_data(filtered_metabolites)
    merged_validation <- validate_metabolite_data(merged_data)
    
    if(!filtered_validation$valid || !merged_validation$valid) {
      stop(paste("Data validation issues found:",
                 paste(c(filtered_validation$messages, 
                         merged_validation$messages), collapse = "\n")))
    }
    
    # Clean and process data
    # Ensure Tissue.type is character
    filtered_metabolites$Tissue.type <- as.character(filtered_metabolites$Tissue.type)
    merged_data$Tissue.type <- as.character(merged_data$Tissue.type)
    
    # Process root and shoot metabolites
    root_metabolites <- filtered_metabolites[filtered_metabolites$Tissue.type == "R", ]
    leaf_metabolites <- filtered_metabolites[filtered_metabolites$Tissue.type == "L", ]
    
    # Identify crosstalk metabolites
    crosstalk_metabolites <- merge(root_metabolites, leaf_metabolites, 
                                   by = "Metabolite", suffixes = c("_root", "_leaf"))
    
    if(nrow(crosstalk_metabolites) == 0) {
      stop("No crosstalk metabolites found")
    }
    
    # Filter merged data
    filtered_merged_data <- merged_data[merged_data$Metabolite %in% 
                                          crosstalk_metabolites$Metabolite, ]
    
    # Ensure Metabolite_Value is numeric
    filtered_merged_data$Metabolite_Value <- as.numeric(filtered_merged_data$Metabolite_Value)
    
    # Add tissue prefix
    filtered_merged_data$Metabolite <- paste0(
      substr(filtered_merged_data$Tissue.type, 1, 1),
      "_",
      filtered_merged_data$Metabolite
    )
    
    # Create pivot data safely
    pivot_data <- filtered_merged_data %>%
      select(Vac_id, Genotype, Entry, Batch, Treatment, Replication, Day, 
             Metabolite, Metabolite_Value) %>%
      pivot_wider(
        id_cols = c(Vac_id, Genotype, Entry, Batch, Treatment, Replication, Day),
        names_from = Metabolite,
        values_from = Metabolite_Value,
        values_fill = list(Metabolite_Value = NA)
      )
    
    # Convert all non-ID columns to numeric
    numeric_cols <- setdiff(colnames(pivot_data), 
                            c("Vac_id", "Genotype", "Entry", "Batch", "Treatment", "Replication", "Day"))
    pivot_data[numeric_cols] <- lapply(pivot_data[numeric_cols], as.numeric)
    
    # Create processing summary
    processing_summary <- list(
      n_root_metabolites = nrow(root_metabolites),
      n_leaf_metabolites = nrow(leaf_metabolites),
      n_crosstalk = nrow(crosstalk_metabolites)
    )
    
    return(list(
      pivot_data = pivot_data,
      crosstalk_metabolites = crosstalk_metabolites,
      processing_summary = processing_summary
    ))
    
  }, error = function(e) {
    stop("Error in data processing: ", e$message)
  })
}

# Execute processing with error catching
cat("Starting data processing...\n")
processed_data <- load_process_data()

# Save processed data
saveRDS(processed_data, file.path(base_output_path, "processed_data.rds"))

# Print summary
cat("Data processing completed.\n",
    "Number of crosstalk metabolites:", 
    nrow(processed_data$crosstalk_metabolites), "\n",
    "Results saved in:", base_output_path, "\n")



# Section 2: Bayesian Network & Bootstrap 
library(bnlearn)
library(parallel)
library(foreach)
library(doParallel)
library(stats)
library(tidyverse)
library(igraph)

# Clean up any existing parallel connections
if(exists("cluster")) {
  stopCluster(cluster)
  rm(cluster)
}
closeAllConnections()

# Setup parallel processing with explicit cleanup
n_cores <- min(detectCores() - 1, 4)  # Limit to max 4 cores for safety
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

# Modified data loading and preprocessing
load_and_prepare_data <- function(base_output_path) {
  tryCatch({
    processed_data <- readRDS(file.path(base_output_path, "processed_data.rds"))
    pivot_data <- processed_data$pivot_data
    
    # Print data dimensions for debugging
    cat("Loaded data dimensions:", dim(pivot_data), "\n")
    
    # Ensure numeric conversion and handle NAs
    data_for_bn <- pivot_data[, -(1:7)] %>%
      mutate(across(everything(), as.numeric)) %>%
      mutate(across(everything(), ~ifelse(is.infinite(.), NA, .)))
    
    # Remove columns with all NAs
    data_for_bn <- data_for_bn[, colSums(is.na(data_for_bn)) < nrow(data_for_bn)]
    cat("After NA handling dimensions:", dim(data_for_bn), "\n")
    
    # Impute remaining NAs with median
    data_for_bn <- as.data.frame(apply(data_for_bn, 2, function(x) {
      ifelse(is.na(x), median(x, na.rm = TRUE), x)
    }))
    
    # Split into root and leaf data
    root_data <- data_for_bn[, grep("^R_", colnames(data_for_bn), value = TRUE)]
    leaf_data <- data_for_bn[, grep("^L_", colnames(data_for_bn), value = TRUE)]
    
    cat("Root data dimensions:", dim(root_data), "\n")
    cat("Leaf data dimensions:", dim(leaf_data), "\n")
    
    if(ncol(root_data) == 0 || ncol(leaf_data) == 0) {
      stop("No valid columns found in either root or leaf data")
    }
    
    return(list(
      root_data = root_data,
      leaf_data = leaf_data,
      crosstalk_metabolites = processed_data$crosstalk_metabolites
    ))
  }, error = function(e) {
    stop("Error in data preparation: ", e$message)
  })
}

# Non-parametric bootstrap analysis
bootstrap_network_analysis <- function(data, n_boots = 5000) {  
  if(ncol(data) == 0 || nrow(data) == 0) {
    stop("Empty dataset provided")
  }
  
  cat("Bootstrap analysis for dataset:", dim(data), "\n")
  
  # Initialize storage
  edge_counts <- matrix(0, ncol(data), ncol(data))
  colnames(edge_counts) <- rownames(edge_counts) <- colnames(data)
  
  # Perform bootstrap iterations
  successful_boots <- 0
  boot_results <- matrix(0, ncol(data), ncol(data))
  
  for(i in 1:n_boots) {
    tryCatch({
      # Non-parametric bootstrap sampling
      boot_indices <- sample(nrow(data), replace = TRUE)
      boot_data <- data[boot_indices, ]
      
      # Learn network structure
      boot_net <- hc(boot_data)
      
      # Update results
      boot_results <- boot_results + amat(boot_net)
      successful_boots <- successful_boots + 1
      
    }, error = function(e) {
      cat("Warning: Bootstrap iteration", i, "failed:", e$message, "\n")
    })
  }
  
  if(successful_boots == 0) {
    stop("All bootstrap iterations failed")
  }
  
  # Calculate edge probabilities and significance
  edge_prob <- boot_results / successful_boots
  edge_pvals <- 1 - edge_prob
  adj_pvals <- p.adjust(as.vector(edge_pvals), method = "BH")
  significant_edges <- matrix(adj_pvals < 0.05, nrow = ncol(data))
  
  cat("Bootstrap analysis completed with", successful_boots, "iterations\n")
  
  return(list(
    edge_prob = edge_prob,
    significant_edges = significant_edges,
    adj_pvals = matrix(adj_pvals, nrow = ncol(data))
  ))
}

# Fixed null model generation
generate_null_networks <- function(data, n_permutations = 5000) {  
  if(ncol(data) == 0 || nrow(data) == 0) {
    stop("Empty dataset provided")
  }
  
  cat("Generating null model for dataset with dimensions:", dim(data), "\n")
  
  # Calculate observed network
  observed_network <- hc(data)
  cat("Observed network nodes:", length(nodes(observed_network)), "\n")
  cat("Observed network arcs:", nrow(arcs(observed_network)), "\n")
  
  # Store results
  observed_metrics <- list(
    n_edges = nrow(arcs(observed_network)),
    avg_markov_blanket = mean(sapply(nodes(observed_network), 
                                     function(n) length(mb(observed_network, n))))
  )
  
  # Generate null metrics
  null_metrics <- matrix(NA, nrow = n_permutations, ncol = length(observed_metrics))
  colnames(null_metrics) <- names(observed_metrics)
  
  cat("Computing null models...\n")
  successful_perms <- 0
  
  for(i in 1:n_permutations) {
    tryCatch({
      # Permute within each column
      null_data <- as.data.frame(apply(data, 2, sample))
      colnames(null_data) <- colnames(data)
      
      # Learn network
      null_network <- hc(null_data)
      
      # Calculate metrics
      null_metrics[i,] <- c(
        n_edges = nrow(arcs(null_network)),
        avg_markov_blanket = mean(sapply(nodes(null_network), 
                                         function(n) length(mb(null_network, n))))
      )
      
      successful_perms <- successful_perms + 1
      if(i %% 5 == 0) cat("Completed", i, "permutations\n")
      
    }, error = function(e) {
      cat("Warning: Permutation", i, "failed:", e$message, "\n")
    })
  }
  
  # Remove failed iterations
  null_metrics <- na.omit(null_metrics)
  
  if(nrow(null_metrics) == 0) {
    stop("All null model iterations failed")
  }
  
  cat("Successfully generated", nrow(null_metrics), "null models\n")
  
  # Calculate empirical p-values
  p_values <- sapply(1:ncol(null_metrics), function(i) {
    mean(null_metrics[,i] >= observed_metrics[[i]], na.rm = TRUE)
  })
  names(p_values) <- names(observed_metrics)
  
  # Calculate robust effect sizes
  effect_sizes <- sapply(1:ncol(null_metrics), function(i) {
    (observed_metrics[[i]] - median(null_metrics[,i], na.rm = TRUE)) / 
      IQR(null_metrics[,i], na.rm = TRUE)
  })
  names(effect_sizes) <- names(observed_metrics)
  
  return(list(
    observed = observed_metrics,
    null_mean = apply(null_metrics, 2, median, na.rm = TRUE),
    null_sd = apply(null_metrics, 2, IQR, na.rm = TRUE),
    p_values = p_values,
    effect_sizes = effect_sizes,
    n_permutations = nrow(null_metrics)
  ))
}

# Main execution block
cat("Starting enhanced Bayesian Network analysis...\n")
base_output_path <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"

tryCatch({
  # Load and prepare data with diagnostics
  cat("Loading and preparing data...\n")
  prepared_data <- load_and_prepare_data(base_output_path)
  
  # Perform bootstrap analysis
  cat("Analyzing root network...\n")
  root_network_results <- bootstrap_network_analysis(prepared_data$root_data)
  
  cat("Analyzing leaf network...\n")
  leaf_network_results <- bootstrap_network_analysis(prepared_data$leaf_data)
  
  # Generate null models
  cat("Generating null models for root network...\n")
  root_null_results <- generate_null_networks(prepared_data$root_data)
  
  cat("Generating null models for leaf network...\n")
  leaf_null_results <- generate_null_networks(prepared_data$leaf_data)
  
  # Combine all results
  network_results <- list(
    root_network = root_network_results,
    leaf_network = leaf_network_results,
    root_null = root_null_results,
    leaf_null = leaf_null_results
  )
  
  # Save results
  saveRDS(network_results, file.path(base_output_path, "network_results.rds"))
  
  # Save summary statistics
  null_model_summary <- data.frame(
    Network = c("Root", "Leaf"),
    Observed_Edges = c(root_null_results$observed$n_edges, 
                       leaf_null_results$observed$n_edges),
    Null_Mean_Edges = c(root_null_results$null_mean["n_edges"], 
                        leaf_null_results$null_mean["n_edges"]),
    Effect_Size = c(root_null_results$effect_sizes["n_edges"], 
                    leaf_null_results$effect_sizes["n_edges"]),
    P_Value = c(root_null_results$p_values["n_edges"], 
                leaf_null_results$p_values["n_edges"])
  )
  
  write.csv(null_model_summary, 
            file.path(base_output_path, "null_model_summary.csv"), 
            row.names = FALSE)
  
  cat("Analysis completed successfully.\n")
  cat("Results saved to:", base_output_path, "\n")
  
}, error = function(e) {
  cat("Error in main execution:", conditionMessage(e), "\n")
  cat("Stack trace:\n")
  print(sys.calls())
}, finally = {
  # Clean up parallel backend
  stopCluster(cluster)
  closeAllConnections()
  gc()
})


# Section 3: Temporal Analysis 
library(dtw)
library(tidyverse)

# Define base path
base_output_path <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"

# Load previously processed data
processed_data <- readRDS(file.path(base_output_path, "processed_data.rds"))
network_results <- readRDS(file.path(base_output_path, "network_results.rds"))


# Enhanced debug function with data validation
check_data_structure <- function(pivot_data) {
  tryCatch({
    cat("Data structure validation:\n")
    cat("Dimensions:", dim(pivot_data), "\n")
    cat("Column names:", paste(head(colnames(pivot_data), 5), collapse=", "), "...\n")
    days <- sort(unique(pivot_data$Day))
    cat("Time points:", paste(days, collapse=", "), "\n")
    cat("Minimum samples per time point:", 
        min(table(pivot_data$Day)), "\n")
    
    # Validate time series structure
    if(length(days) < 2) {
      stop("Insufficient time points for temporal analysis")
    }
    return(TRUE)
  }, error = function(e) {
    cat("Data structure validation failed:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Enhanced temporal coherence calculation with FDR correction
calculate_temporal_coherence <- function(pivot_data, n_permutations = 20) {
  tryCatch({
    # Validate data structure
    if(!check_data_structure(pivot_data)) {
      stop("Invalid data structure for temporal analysis")
    }
    
    # Separate root and leaf data
    root_cols <- grep("^R_", colnames(pivot_data), value = TRUE)
    leaf_cols <- grep("^L_", colnames(pivot_data), value = TRUE)
    
    # Get metabolite names without prefixes
    root_mets <- gsub("^R_", "", root_cols)
    leaf_mets <- gsub("^L_", "", leaf_cols)
    common_metabolites <- intersect(root_mets, leaf_mets)
    
    cat("Found", length(root_cols), "root metabolites\n")
    cat("Found", length(leaf_cols), "leaf metabolites\n")
    cat("Found", length(common_metabolites), "common metabolites\n")
    
    # Initialize results storage
    all_results <- list()
    
    # Process each metabolite
    for(met in common_metabolites) {
      cat("Processing metabolite:", met, "\n")
      
      root_name <- paste0("R_", met)
      leaf_name <- paste0("L_", met)
      
      # Extract and validate time series
      leaf_ts <- tapply(pivot_data[[leaf_name]], pivot_data$Day, 
                        median, na.rm = TRUE)  # Using median for robustness
      root_ts <- tapply(pivot_data[[root_name]], pivot_data$Day, 
                        median, na.rm = TRUE)  # Using median for robustness
      
      # Debug information
      cat("Time series lengths - Leaf:", length(leaf_ts), 
          "Root:", length(root_ts), "\n")
      cat("NA check - Leaf:", sum(is.na(leaf_ts)), 
          "Root:", sum(is.na(root_ts)), "\n")
      
      # Skip if insufficient data
      if(length(leaf_ts) < 2 || length(root_ts) < 2 || 
         any(is.na(leaf_ts)) || any(is.na(root_ts)) ||
         length(leaf_ts) != length(root_ts)) {
        cat("Skipping", met, "due to insufficient data\n")
        next
      }
      
      # Calculate DTW distance
      actual_dtw <- dtw(leaf_ts, root_ts, distance.only = TRUE)$distance
      
      # Permutation test
      perm_distances <- numeric(n_permutations)
      for(i in seq_along(perm_distances)) {
        perm_root <- sample(root_ts)
        perm_distances[i] <- dtw(leaf_ts, perm_root, 
                                 distance.only = TRUE)$distance
      }
      
      # Calculate robust statistics
      p_value <- mean(perm_distances <= actual_dtw)
      effect_size <- (median(perm_distances) - actual_dtw) / 
        IQR(perm_distances)  # Robust effect size
      
      # Calculate rank correlation
      temp_cor <- cor(leaf_ts, root_ts, 
                      method = "spearman",    # Non-parametric correlation
                      use = "pairwise.complete.obs")
      
      # Store results
      all_results[[met]] <- data.frame(
        Metabolite = met,
        DTW_Distance = actual_dtw,
        Temporal_Correlation = temp_cor,
        P_Value = p_value,
        Effect_Size = effect_size,
        N_Timepoints = length(leaf_ts)
      )
    }
    
    # Combine results
    if(length(all_results) == 0) {
      stop("No valid results obtained for any metabolite")
    }
    
    coherence_results <- do.call(rbind, all_results)
    
    # Add FDR-corrected p-values
    coherence_results$Adjusted_P_Value <- p.adjust(coherence_results$P_Value, 
                                                   method = "BH")
    
    # Add significance levels
    coherence_results$Significance <- cut(coherence_results$Adjusted_P_Value,
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "ns"))
    
    return(coherence_results)
  }, error = function(e) {
    cat("Error details:", conditionMessage(e), "\n")
    stop("Error in temporal coherence calculation: ", e$message)
  })
}

# Enhanced analysis function with network integration
analyze_temporal_patterns <- function(coherence_results, network_results) {
  # Integrate temporal coherence with network structure
  temporal_network_overlap <- data.frame(
    coherence_results,
    Root_Node = paste0("R_", coherence_results$Metabolite),
    Leaf_Node = paste0("L_", coherence_results$Metabolite)
  )
  
  # Add network metrics
  temporal_network_overlap$Root_Degree <- sapply(temporal_network_overlap$Root_Node,
                                                 function(node) sum(network_results$root_network$edge_prob[node,] > 0.5))
  temporal_network_overlap$Leaf_Degree <- sapply(temporal_network_overlap$Leaf_Node,
                                                 function(node) sum(network_results$leaf_network$edge_prob[node,] > 0.5))
  
  return(temporal_network_overlap)
}

# Execute temporal analysis with enhanced error handling
cat("Starting temporal analysis...\n")

tryCatch({
  # Calculate temporal coherence
  temporal_coherence <- calculate_temporal_coherence(processed_data$pivot_data)
  
  # Analyze patterns
  temporal_patterns <- analyze_temporal_patterns(temporal_coherence, network_results)
  
  # Create comprehensive summary
  temporal_summary <- list(
    n_significant_temporal = sum(temporal_patterns$Adjusted_P_Value < 0.05),
    mean_effect_size = median(temporal_patterns$Effect_Size, na.rm = TRUE),
    temporal_network_overlap = sum(temporal_patterns$Root_Degree > 0 & 
                                     temporal_patterns$Leaf_Degree > 0),
    temporal_coherence = temporal_coherence,
    temporal_patterns = temporal_patterns
  )
  
  # Save results
  saveRDS(temporal_summary, file.path(base_output_path, "temporal_analysis.rds"))
  
  # Save detailed results for publication
  write.csv(temporal_patterns, 
            file.path(base_output_path, "temporal_patterns_detailed.csv"),
            row.names = FALSE)
  
  # Print summary
  cat("\nTemporal Analysis Summary:\n")
  cat("Significant temporal associations:", temporal_summary$n_significant_temporal, "\n")
  cat("Median effect size:", round(temporal_summary$mean_effect_size, 3), "\n")
  cat("Network-temporal overlaps:", temporal_summary$temporal_network_overlap, "\n")
  cat("Results saved to:", file.path(base_output_path, "temporal_analysis.rds"), "\n")
  
}, error = function(e) {
  cat("Error in main execution:", conditionMessage(e), "\n")
  cat("Please check the data structure and requirements.\n")
})

# Explicit garbage collection
gc()



# Section 4: Network Metrics & Validation 
library(tidyverse)
library(igraph)
library(reshape2)
library(RColorBrewer)
library(coin)     # For permutation tests
library(boot)     # For bootstrapping
library(Hmisc)    # For robust correlations

# Define file paths
BASE_DIR <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_new_crosstalk5_V5_5000"
INPUT_FILE <- file.path(BASE_DIR, "network_results.rds")
OUTPUT_DIR <- file.path(BASE_DIR, "network_analysis_results")

# Debug message function
debug_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

# Calculate Network Metrics with Statistical Validation
calculate_network_metrics <- function(network_data, type = c("root", "leaf")) {
  type <- match.arg(type)
  
  tryCatch({
    cat(sprintf("\n=== Processing %s network ===\n", type))
    cat("\nCreating adjacency matrix...\n")
    
    if(is.null(network_data$edge_prob)) {
      stop("No edge probability data found")
    }
    
    adj_matrix <- network_data$edge_prob > 0.5
    cat(sprintf("Adjacency matrix dimensions: %d x %d\n", nrow(adj_matrix), ncol(adj_matrix)))
    
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "max", weighted = TRUE)
    n_edges <- gsize(g)  # Updated from ecount
    cat(sprintf("Graph created with %d vertices and %d edges\n", gorder(g), n_edges))  # Updated from vcount
    
    degrees <- degree(g, mode = "all")
    mean_deg <- mean(degrees)
    median_deg <- median(degrees)
    
    g_undir <- as_undirected(g, mode = "collapse")  # Updated from as.undirected
    comm <- cluster_louvain(g_undir)
    
    metrics <- list(
      network_type = type,
      n_nodes = gorder(g),  # Updated from vcount
      n_edges = n_edges,
      density = edge_density(g),
      mean_degree = mean_deg,
      median_degree = median_deg,
      transitivity = transitivity(g, type = "global"),
      modularity = modularity(comm),
      n_components = components(g)$no,
      mean_path = mean_distance(g, directed = FALSE),
      diameter = diameter(g, directed = FALSE),
      assortativity = assortativity_degree(g, directed = FALSE)
    )
    
    metrics$graph <- g
    
    return(metrics)
    
  }, error = function(e) {
    cat(sprintf("\nERROR in %s network analysis: %s\n", type, e$message))
    return(NULL)
  })
}



# Compare networks with improved data handling
compare_networks <- function(root_metrics, leaf_metrics) {
  cat("\n=== Comparing Networks ===\n")
  
  if(is.null(root_metrics) || is.null(leaf_metrics)) {
    cat("ERROR: Missing metrics, cannot compare\n")
    return(NULL)
  }
  
  # Get common numeric metrics (excluding network_type and graph object)
  root_metrics_clean <- root_metrics[!names(root_metrics) %in% c("network_type", "graph")]
  leaf_metrics_clean <- leaf_metrics[!names(leaf_metrics) %in% c("network_type", "graph")]
  
  # Create comparison dataframe
  comparison <- data.frame(
    Metric = names(root_metrics_clean),
    Root = unlist(root_metrics_clean),
    Leaf = unlist(leaf_metrics_clean)
  ) %>%
    mutate(
      Difference = Root - Leaf,
      Ratio = Root / Leaf,
      PercentDiff = ((Root - Leaf) / ((Root + Leaf)/2)) * 100
    )
  
  return(comparison)
}


# Calculate hub metrics for tissue-specific analysis
calculate_hub_metrics <- function(g, type) {
  tryCatch({
    betweenness <- betweenness(g, normalized = TRUE)
    closeness <- closeness(g, normalized = TRUE)
    eigen <- eigen_centrality(g)$vector
    
    n_hubs <- ceiling(gorder(g) * 0.1)  # Updated from vcount
    hub_indices <- order(betweenness, decreasing = TRUE)[1:n_hubs]
    
    hub_metrics <- data.frame(
      Tissue = type,
      Node = V(g)$name[hub_indices],
      Betweenness = betweenness[hub_indices],
      Closeness = closeness[hub_indices],
      Eigenvector = eigen[hub_indices]
    )
    
    return(hub_metrics)
  }, error = function(e) {
    cat(sprintf("\nERROR in hub metrics calculation: %s\n", e$message))
    return(NULL)
  })
}


calculate_module_metrics <- function(g, type) {
  tryCatch({
    # Convert to undirected if needed
    if(is_directed(g)) {  # Updated from is.directed()
      g <- as_undirected(g, mode="collapse")
    }
    
    # Ensure graph is connected and convert sizes to numeric
    components <- components(g)
    if(components$no > 1) {
      largest_comp <- which.max(components$csize)
      g <- induced_subgraph(g, which(components$membership == largest_comp))
    }
    
    # Detect communities using multiple algorithms for robustness
    communities <- tryCatch({
      cluster_louvain(g)
    }, error = function(e) {
      tryCatch({
        cluster_fast_greedy(g)
      }, error = function(e) {
        cluster_walktrap(g)
      })
    })
    
    if(length(communities) == 0) {
      stop("No communities detected")
    }
    
    # Ensure sizes are numeric
    community_sizes <- as.numeric(sizes(communities))
    
    # Calculate metrics for each module
    module_metrics <- data.frame(
      Tissue = type,
      Module = seq_len(length(communities)),
      Size = community_sizes  # Using explicitly numeric sizes
    )
    # Rest of the function remains the same...
    
    # Add summary statistics with explicit numeric conversion
    cat(sprintf("\nModule analysis for %s tissue:", type))
    cat(sprintf("\nNumber of modules: %d", nrow(module_metrics)))
    cat(sprintf("\nMean module size: %.2f", mean(as.numeric(module_metrics$Size), na.rm=TRUE)))
    cat(sprintf("\nMean internal density: %.3f", mean(module_metrics$Internal_Density, na.rm=TRUE)))
    cat("\n")
    
    return(module_metrics)
    
  }, error = function(e) {
    cat(sprintf("\nERROR in module metrics calculation for %s: %s\n", type, e$message))
    return(data.frame(
      Tissue = type,
      Module = 1,
      Size = as.numeric(gorder(g)),
      Internal_Density = edge_density(g),
      Isolation = 1,
      Avg_Path_Length = mean_distance(g),
      Transitivity = transitivity(g)
    ))
  })
}


# Calculate integration metrics between tissues
calculate_integration_metrics <- function(root_g, leaf_g) {
  tryCatch({
    integration_metrics <- data.frame(
      Metric = c("Network_Assortativity", "Cross_Tissue_Density", "Path_Length_Ratio"),
      Root_Value = c(
        assortativity_degree(root_g),
        edge_density(root_g),
        mean_distance(root_g)
      ),
      Leaf_Value = c(
        assortativity_degree(leaf_g),
        edge_density(leaf_g),
        mean_distance(leaf_g)
      )
    )
    
    integration_metrics$Integration_Score <- with(integration_metrics,
                                                  ifelse(Root_Value != 0 & Leaf_Value != 0,
                                                         abs(Root_Value - Leaf_Value) / ((Root_Value + Leaf_Value)/2),
                                                         NA))
    
    return(integration_metrics)
  }, error = function(e) {
    cat(sprintf("\nERROR in integration metrics calculation: %s\n", e$message))
    return(NULL)
  })
}

# Generate detailed report with proper file handling
generate_detailed_report <- function(root_metrics, leaf_metrics, comparison, output_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  temp_files <- file.path(OUTPUT_DIR, c("temp_metadata.csv", "temp_root.csv", "temp_leaf.csv", "temp_comparison.csv"))
  
  metadata_df <- data.frame(
    Category = "Analysis Information",
    Metric = c("Analysis Date", "Number of Metrics", "Total Comparisons"),
    Value = c(timestamp, 
              length(root_metrics) - 1,  # Subtract 1 for graph object
              nrow(comparison))
  )
  write.csv(metadata_df, temp_files[1], row.names = FALSE)
  
  root_df <- data.frame(
    Category = "Root Network",
    Metric = names(root_metrics)[names(root_metrics) != "graph"],
    Value = unlist(root_metrics[names(root_metrics) != "graph"])
  )
  write.csv(root_df, temp_files[2], row.names = FALSE)
  
  leaf_df <- data.frame(
    Category = "Leaf Network",
    Metric = names(leaf_metrics)[names(leaf_metrics) != "graph"],
    Value = unlist(leaf_metrics[names(leaf_metrics) != "graph"])
  )
  write.csv(leaf_df, temp_files[3], row.names = FALSE)
  
  write.csv(comparison, temp_files[4], row.names = FALSE)
  
  file.create(output_file)
  for(temp_file in temp_files) {
    file.append(output_file, temp_file)
    file.remove(temp_file)
  }
}

# Main execution function
main <- function() {
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n=== Starting Network Analysis ===\n")
  cat("Input file:", INPUT_FILE, "\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
  
  tryCatch({
    if (!file.exists(INPUT_FILE)) {
      stop(sprintf("Input file not found: %s", INPUT_FILE))
    }
    
    network_results <- readRDS(INPUT_FILE)
    cat("\nNetwork results loaded successfully\n")
    
    root_metrics <- calculate_network_metrics(network_results$root_network, "root")
    leaf_metrics <- calculate_network_metrics(network_results$leaf_network, "leaf")
    
    if (!is.null(root_metrics) && !is.null(leaf_metrics)) {
      cat("\nSaving results...\n")
      
      # Save basic network comparison
      comparison <- compare_networks(root_metrics, leaf_metrics)
      if(!is.null(comparison)) {
        comparison_rounded <- comparison %>%
          mutate(across(where(is.numeric), ~round(., 4)))
        
        write.csv(comparison_rounded, 
                  file.path(OUTPUT_DIR, "network_comparison.csv"), 
                  row.names = FALSE)
        
        # Generate and save additional metrics
        hub_metrics_root <- calculate_hub_metrics(root_metrics$graph, "Root")
        hub_metrics_leaf <- calculate_hub_metrics(leaf_metrics$graph, "Leaf")
        if(!is.null(hub_metrics_root) && !is.null(hub_metrics_leaf)) {
          write.csv(rbind(hub_metrics_root, hub_metrics_leaf),
                    file.path(OUTPUT_DIR, "tissue_hub_metrics.csv"),
                    row.names = FALSE)
        }
        
        module_metrics_root <- calculate_module_metrics(root_metrics$graph, "Root")
        module_metrics_leaf <- calculate_module_metrics(leaf_metrics$graph, "Leaf")
        if(!is.null(module_metrics_root) && !is.null(module_metrics_leaf)) {
          write.csv(rbind(module_metrics_root, module_metrics_leaf),
                    file.path(OUTPUT_DIR, "tissue_module_metrics.csv"),
                    row.names = FALSE)
        }
        
        integration_metrics <- calculate_integration_metrics(root_metrics$graph, leaf_metrics$graph)
        if(!is.null(integration_metrics)) {
          write.csv(integration_metrics,
                    file.path(OUTPUT_DIR, "tissue_integration_metrics.csv"),
                    row.names = FALSE)
        }
        
        # Generate detailed report
        generate_detailed_report(
          root_metrics,
          leaf_metrics,
          comparison_rounded,
          file.path(OUTPUT_DIR, "detailed_network_analysis.csv")
        )
        
        # Save raw results
        saveRDS(list(
          root_metrics = root_metrics,
          leaf_metrics = leaf_metrics,
          comparison = comparison,
          hub_metrics = list(root = hub_metrics_root, leaf = hub_metrics_leaf),
          module_metrics = list(root = module_metrics_root, leaf = module_metrics_leaf),
          integration_metrics = integration_metrics,
          timestamp = Sys.time()
        ), file.path(OUTPUT_DIR, "raw_network_analysis.rds"))
        
        cat("\nAnalysis complete. Results saved in:", OUTPUT_DIR, "\n")
        cat("Main results files:\n")
        cat("1. detailed_network_analysis.csv\n")
        cat("2. tissue_hub_metrics.csv\n")
        cat("3. tissue_module_metrics.csv\n")
        cat("4. tissue_integration_metrics.csv\n")
      }
    }
    
  }, error = function(e) {
    cat("\nERROR in main execution:", conditionMessage(e), "\n")
    cat("Stack trace:\n")
    print(sys.calls())
  })
}

# Execute analysis
main()