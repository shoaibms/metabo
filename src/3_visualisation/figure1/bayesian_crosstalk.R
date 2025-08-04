# -----------------------------------------------------------------------------
# Bayesian Network Analysis of Root-Shoot Metabolite Crosstalk
#
# Description:
# This script performs a Bayesian network analysis to model the interactions
# between metabolites in root and leaf tissues. It identifies crosstalk
# metabolites present in both tissues, constructs separate Bayesian networks
# for each tissue, and then combines them into a single network. The script
# visualizes the combined network in both 2D and 3D, highlighting the
# connections within and between tissues.
#
# The final outputs are:
# - A 2D plot of the interaction network (`root_shoot_interaction_network_2d.png`)
# - An interactive 3D plot of the network (`root_shoot_interaction_network_3d.html`)
#
# -----------------------------------------------------------------------------

# Load required packages
suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(plotly)
  library(reshape2)
  library(bnlearn)
  library(htmlwidgets)
})

# --- Configuration ---

# Hardcoded file paths (as per user request)
filtered_metabolites_file <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/group/merge2/Merged_VIP_hub_name_only.csv"
merged_data_file <- "C:/Users/ms/Desktop/data_chem_3_10/output/results/group/merge2/Merged_VIP_hub.csv"
output_dir <- "C:/Users/ms/Desktop/r/chem_data/final/baysian_root_shoot_plot2"

# Define output file paths
output_2d_plot <- file.path(output_dir, "root_shoot_interaction_network_2d.png")
output_3d_plot <- file.path(output_dir, "root_shoot_interaction_network_3d.html")

# Define constants for analysis
ROOT_PREFIX <- "R_"
LEAF_PREFIX <- "L_"
ROOT_COLOR <- "#1e597d"  # Blue
LEAF_COLOR <- "#1b941b"  # Green
CROSSTALK_COLOR <- "orange"

# Set seed for reproducibility of layouts
set.seed(42)

# --- Utility Functions ---

# Function to check if a directory exists and create it if it doesn't
create_dir_safe <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

# --- Main Script ---

# Create the output directory if it doesn't exist
create_dir_safe(output_dir)

# Validate that input files exist before proceeding
if (!file.exists(filtered_metabolites_file) || !file.exists(merged_data_file)) {
  stop("Input files not found. Please check the specified file paths.")
}

# --- Data Loading and Preprocessing ---

# Load data from CSV files
filtered_metabolites <- read.csv(filtered_metabolites_file)
merged_data <- read.csv(merged_data_file)

# Filter metabolites based on tissue type
root_metabolites <- filtered_metabolites[filtered_metabolites$Tissue.type == "R", ]
leaf_metabolites <- filtered_metabolites[filtered_metabolites$Tissue.type == "L", ]

# Identify crosstalk metabolites that are present in both root and leaf tissues
crosstalk_metabolites <- merge(root_metabolites, leaf_metabolites, by = "Metabolite", suffixes = c("_root", "_leaf"))

# Filter the main dataset for crosstalk metabolites and add tissue prefixes to names
filtered_merged_data <- merged_data[merged_data$Metabolite %in% crosstalk_metabolites$Metabolite, ]
filtered_merged_data$Metabolite <- paste0(substr(filtered_merged_data$Tissue.type, 1, 1), "_", filtered_merged_data$Metabolite)

# Pivot data from long to wide format for network analysis
pivot_data <- dcast(filtered_merged_data, 
                    Vac_id + Genotype + Entry + Batch + Treatment + Replication + Day ~ Metabolite, 
                    value.var = "Metabolite_Value")

# Prepare data for Bayesian network construction by removing metadata columns
metadata_cols <- c("Vac_id", "Genotype", "Entry", "Batch", "Treatment", "Replication", "Day")
data_for_bn <- pivot_data[, !(colnames(pivot_data) %in% metadata_cols)]

# Separate data into root and leaf datasets
root_data <- data_for_bn[, grepl(paste0("^", ROOT_PREFIX), colnames(data_for_bn))]
leaf_data <- data_for_bn[, grepl(paste0("^", LEAF_PREFIX), colnames(data_for_bn))]


# --- Bayesian Network Construction ---

# Learn the network structure for each tissue using the hill-climbing algorithm
bn_structure_root <- hc(root_data)
bn_structure_leaf <- hc(leaf_data)

# Extract arcs (directed edges) from the learned structures
root_arcs <- arcs(bn_structure_root)
leaf_arcs <- arcs(bn_structure_leaf)

# Stop execution if no network structure could be learned
if (nrow(root_arcs) == 0 || nrow(leaf_arcs) == 0) {
  stop("No edges found in the learned network structure. Cannot proceed.")
}

# Create igraph objects from the network structures
graph_root <- graph_from_data_frame(root_arcs, directed = TRUE)
graph_leaf <- graph_from_data_frame(leaf_arcs, directed = TRUE)

# Combine the root and leaf graphs into a single network
all_vertices <- unique(c(V(graph_root)$name, V(graph_leaf)$name))
combined_graph <- make_empty_graph(n = length(all_vertices), directed = TRUE)
V(combined_graph)$name <- all_vertices

# Add edges from both tissue-specific graphs
all_edges <- rbind(as_edgelist(graph_root), as_edgelist(graph_leaf))
combined_graph <- add_edges(combined_graph, edges = as.vector(t(all_edges)))

# Identify and add cross-tissue connections for metabolites common to both tissues
common_metabolites <- intersect(gsub(paste0("^", ROOT_PREFIX), "", V(graph_root)$name), 
                                gsub(paste0("^", LEAF_PREFIX), "", V(graph_leaf)$name))
cross_edges <- data.frame(
  from = paste0(ROOT_PREFIX, common_metabolites),
  to = paste0(LEAF_PREFIX, common_metabolites)
)
combined_graph <- add_edges(combined_graph, edges = as.vector(t(cross_edges)))

# --- Node and Edge Attribute Assignment ---

# Calculate node degrees (a measure of centrality) to determine node sizes
node_degrees <- degree(combined_graph, mode = "all")

# Normalize node degrees to a consistent range for visual scaling
if (length(unique(node_degrees)) > 1) {
  node_scores_norm <- 1 + 9 * (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees))
} else {
  node_scores_norm <- rep(5, length(node_degrees)) # Assign a default size if all degrees are the same
}

# Assign node sizes, making nodes involved in crosstalk larger for emphasis
V(combined_graph)$size <- ifelse(V(combined_graph)$name %in% unlist(cross_edges), 
                                 node_scores_norm * 5,
                                 node_scores_norm)

# Assign node colors based on tissue type
V(combined_graph)$color <- ifelse(grepl(paste0("^", ROOT_PREFIX), V(combined_graph)$name), ROOT_COLOR, LEAF_COLOR)

# Assign edge colors based on the type of connection (within or between tissues)
edge_list <- as_edgelist(combined_graph)
E(combined_graph)$color <- sapply(1:nrow(edge_list), function(i) {
  from_node <- edge_list[i, 1]
  to_node <- edge_list[i, 2]
  is_from_root <- grepl(paste0("^", ROOT_PREFIX), from_node)
  is_to_root <- grepl(paste0("^", ROOT_PREFIX), to_node)
  is_from_leaf <- grepl(paste0("^", LEAF_PREFIX), from_node)
  is_to_leaf <- grepl(paste0("^", LEAF_PREFIX), to_node)
  
  if (is_from_root && is_to_root) {
    return(ROOT_COLOR)      # Root to Root
  } else if (is_from_leaf && is_to_leaf) {
    return(LEAF_COLOR)      # Leaf to Leaf
  } else {
    return(CROSSTALK_COLOR) # Cross-tissue connections
  }
})

# --- 2D Network Visualization ---

tryCatch({
  g_2d <- ggraph(combined_graph, layout = "fr") +
    geom_edge_link(aes(color = I(color)),
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(2, 'mm'),
                   alpha = 0.7,
                   edge_width = 0.5) +
    geom_node_point(aes(size = size, color = I(color)), alpha = 0.7) +
    scale_size_continuous(range = c(1, 20)) +
    theme_void() +
    labs(title = "Root-Shoot Interaction Network - Crosstalk Molecular Features") +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          legend.position = "none") +
    annotate("text", x = Inf, y = -Inf,
             label = paste("Nodes: Blue = Root, Green = Leaf",
                           "Edges: Blue = Root-Root, Green = Leaf-Leaf, Orange = Root-Leaf",
                           sep = "\n"),
             hjust = 1.1, vjust = -1, size = 5, color = "black")
  
  ggsave(output_2d_plot, g_2d, width = 15, height = 15, dpi = 300)
}, error = function(e) {
  message("Error in 2D plot creation: ", e$message)
})


# --- 3D Network Visualization ---

tryCatch({
  layout_3d <- layout_with_fr(combined_graph, dim = 3)
  
  node_data <- data.frame(
    x = layout_3d[,1],
    y = layout_3d[,2],
    z = layout_3d[,3],
    name = V(combined_graph)$name,
    size = V(combined_graph)$size,
    type = ifelse(grepl(paste0("^", ROOT_PREFIX), V(combined_graph)$name), "Root", "Leaf")
  )
  node_data$color <- ifelse(node_data$type == "Root", ROOT_COLOR, LEAF_COLOR)
  
  edge_list_3d <- as_edgelist(combined_graph)
  edge_segments <- lapply(1:nrow(edge_list_3d), function(i) {
    from_idx <- which(V(combined_graph)$name == edge_list_3d[i, 1])
    to_idx <- which(V(combined_graph)$name == edge_list_3d[i, 2])
    
    from_node <- V(combined_graph)$name[from_idx]
    to_node <- V(combined_graph)$name[to_idx]
    
    edge_type <- if (grepl(paste0("^", ROOT_PREFIX), from_node) && grepl(paste0("^", ROOT_PREFIX), to_node)) {
      "Root-Root"
    } else if (grepl(paste0("^", LEAF_PREFIX), from_node) && grepl(paste0("^", LEAF_PREFIX), to_node)) {
      "Leaf-Leaf"
    } else {
      "Root-Leaf"
    }
    
    data.frame(
      x = c(layout_3d[from_idx, 1], layout_3d[to_idx, 1], NA),
      y = c(layout_3d[from_idx, 2], layout_3d[to_idx, 2], NA),
      z = c(layout_3d[from_idx, 3], layout_3d[to_idx, 3], NA),
      type = edge_type
    )
  })
  edge_data <- do.call(rbind, edge_segments)
  
  p_3d <- plot_ly() %>%
    add_trace(
      data = edge_data,
      x = ~x, y = ~y, z = ~z,
      type = 'scatter3d', mode = 'lines',
      line = list(width = 1),
      color = ~type,
      colors = c("Root-Root" = ROOT_COLOR, "Leaf-Leaf" = LEAF_COLOR, "Root-Leaf" = CROSSTALK_COLOR),
      name = ~type,
      legendgroup = ~type,
      hoverinfo = 'none'
    ) %>%
    add_trace(
      data = node_data,
      x = ~x, y = ~y, z = ~z,
      type = 'scatter3d', mode = 'markers',
      marker = list(size = ~size, color = ~color),
      text = ~name,
      hoverinfo = 'text',
      name = ~type,
      legendgroup = ~type
    ) %>%
    layout(
      scene = list(
        xaxis = list(title = '', showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
        yaxis = list(title = '', showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
        zaxis = list(title = '', showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE)
      ),
      title = "3D Root-Shoot Interaction Network",
      legend = list(title = list(text = 'Node and Edge Types'))
    )
  
  saveWidget(p_3d, output_3d_plot)
}, error = function(e) {
  message("Error in 3D plot creation: ", e$message)
})

# --- Final Output ---

print(paste("2D network visualization saved to:", output_2d_plot))
print(paste("3D network visualization saved to:", output_3d_plot))
