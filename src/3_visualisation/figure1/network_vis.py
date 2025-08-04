import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

# File paths
leaf_data_path = r'C:\Users\ms\Desktop\data_chem_3_10\data\data\n_p_l.csv'
root_data_path = r'C:\Users\ms\Desktop\data_chem_3_10\data\data\n_p_r.csv'
vip_file_path = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\spearman\VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv'
output_dir = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\spearman\filtered_networks_full_PLOT'

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Load VIP file
vip_df = pd.read_csv(vip_file_path)
vip_metabolites = vip_df['Metabolite'].tolist()

def filter_treatment_and_plot(data_path, output_file, title, tissue_color):
    """
    Filters Treatment 1 data, generates a network plot using the full dataset, 
    and highlights VIP metabolites with node sizes adjusted based on node degrees.
    """
    # Load data
    data = pd.read_csv(data_path)
    
    # Filter for Treatment 1
    data = data[data['Treatment'] == 1]
    
    # Extract metabolite columns (assuming these start at column 11)
    data = data[data.columns[10:]]
    
    # Calculate correlation matrix
    corr_matrix = data.corr(method='spearman')
    corr_matrix = corr_matrix[(corr_matrix.abs() > 0.7) & (corr_matrix != 1)]  # Filter correlations

    # Build network
    G = nx.Graph()
    for metabolite1 in corr_matrix.columns:
        for metabolite2 in corr_matrix.index:
            correlation_value = corr_matrix.loc[metabolite1, metabolite2]
            if not pd.isnull(correlation_value):
                G.add_edge(metabolite1, metabolite2, weight=correlation_value)

    # Calculate degrees
    degrees = dict(G.degree())

    # Option 1: Logarithmic Scaling
    degrees_log = {node: np.log(degree + 1) for node, degree in degrees.items()}
    
    # Scaling factors (adjust as needed)
    non_vip_scale = 10
    vip_scale = 15

    # Apply scaling to node sizes
    non_vip_node_sizes = [degrees_log[node] * non_vip_scale for node in G.nodes if node not in vip_metabolites]
    vip_node_sizes = [degrees_log[node] * vip_scale for node in G.nodes if node in vip_metabolites]

    # Separate VIP and non-VIP nodes
    vip_nodes = [node for node in G.nodes if node in vip_metabolites]
    non_vip_nodes = [node for node in G.nodes if node not in vip_metabolites]

    # Visualize network
    pos = nx.spring_layout(G, k=0.3, iterations=30, seed=42)############################5000
    plt.figure(figsize=(6, 6))

    # Draw nodes with sizes based on degrees
    nx.draw_networkx_nodes(
        G, pos, nodelist=non_vip_nodes, node_color=tissue_color,
        node_size=non_vip_node_sizes, alpha=0.6
    )
    nx.draw_networkx_nodes(
        G, pos, nodelist=vip_nodes, node_color="#2df7b7",
        node_size=vip_node_sizes, alpha=0.50
    )

    # Draw edges
    positive_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] > 0]
    negative_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] < 0]

    nx.draw_networkx_edges(G, pos, edgelist=positive_edges, edge_color="#a69050", width=0.15, alpha=0.6)
    nx.draw_networkx_edges(G, pos, edgelist=negative_edges, edge_color="#5d8ba3", width=0.15, alpha=0.6)

    # Add legend
    legend_elements = [
        mpatches.Patch(color=tissue_color, label='Molecular Features'),
        mpatches.Patch(color="#2df7b7", label='VIP Features'),
        mlines.Line2D([], [], color="#a69050", linewidth=2, label='Positive Correlation'),
        mlines.Line2D([], [], color="#5d8ba3", linewidth=2, label='Negative Correlation'),
        mlines.Line2D([], [], color="black", marker='o', linestyle='None',
                      markersize=10, label='Node Size ~ Degree')
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize='medium')
    plt.title(title, fontsize=16)
    plt.axis('off')
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

# Generate plots for root and leaf
filter_treatment_and_plot(
    leaf_data_path,
    os.path.join(output_dir, 'leaf_network_full_treatment1.png'),
    'Leaf Network',
    '#2ecc71'
)
filter_treatment_and_plot(
    root_data_path,
    os.path.join(output_dir, 'root_network_full_treatment1.png'),
    'Root Network',
    '#33d6d3'
)

print("Plots generated and saved.")
