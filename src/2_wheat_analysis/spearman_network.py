# Section 0: Import necessary libraries
import joblib
import pandas as pd
import numpy as np
import os
import networkx as nx
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from openpyxl import Workbook
from datetime import datetime
import matplotlib.lines as mlines
from matplotlib import cm
from scipy.stats import spearmanr
import community as community_louvain
import plotly.graph_objs as go 
from statistics import median
import multiprocessing
from functools import partial
import pickle
import warnings
warnings.filterwarnings('ignore')

# Paths to the data files and output directory
leaf_data_path = r'C:\Users\USER\Desktop\data_chem_3_10\data\data\n_p_l.csv'
root_data_path = r'C:\Users\USER\Desktop\data_chem_3_10\data\data\n_p_r.csv'
output_dir = r'C:\Users\USER\Desktop\data_chem_3_10\output\results\spearman\network2_plot4E'
vip_file_path = r'C:\Users\USER\Desktop\data_chem_3_10\output\results\spearman\VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv'

# Load the VIP file (Ensure this is defined before using it)
vip_file_df = pd.read_csv(vip_file_path)

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

def calculate_spearman_correlations(data):
    """Calculate Spearman correlations with error handling"""
    try:
        correlation_matrix, _ = spearmanr(data)
        correlation_matrix_df = pd.DataFrame(correlation_matrix, index=data.columns, columns=data.columns)
        return correlation_matrix_df
    except Exception as e:
        print(f"Error in correlation calculation: {e}")
        return None

def filter_correlations(corr_matrix, threshold=0.7):
    """Filter correlations by threshold with error handling"""
    try:
        strong_corr = corr_matrix[(corr_matrix.abs() > threshold) & (corr_matrix != 1)]
        return strong_corr
    except Exception as e:
        print(f"Error in correlation filtering: {e}")
        return None

def construct_network_with_sign(corr_matrix):
    """Construct network preserving correlation signs"""
    G = nx.Graph()
    for metabolite1 in corr_matrix.columns:
        for metabolite2 in corr_matrix.index:
            correlation_value = corr_matrix.loc[metabolite1, metabolite2]
            if not pd.isnull(correlation_value):
                G.add_edge(metabolite1, metabolite2, weight=correlation_value)
    return G

def calculate_extended_network_metrics(G):
    """Calculate comprehensive network metrics"""
    metrics = {}
    
    # Basic network properties
    metrics['nodes'] = nx.number_of_nodes(G)
    metrics['edges'] = nx.number_of_edges(G)
    metrics['density'] = nx.density(G)
    
    # Clustering and community metrics
    metrics['transitivity'] = nx.transitivity(G)
    
    # Calculate modularity using community detection
    G_abs = nx.Graph()
    for u, v, data in G.edges(data=True):
        G_abs.add_edge(u, v, weight=abs(data.get('weight', 1.0)))
    try:
        communities = community_louvain.best_partition(G_abs)
        metrics['modularity'] = community_louvain.modularity(communities, G_abs)
    except ValueError as e:
        print(f"Warning: Modularity calculation failed with error: {e}")
        metrics['modularity'] = 0.0
    
    # Degree metrics
    degrees = [d for n, d in G.degree()]
    metrics['median_degree'] = median(degrees)
    
    # Eigenvector centrality
    try:
        eig_cent = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
        metrics['eigenvector_centrality_mean'] = np.mean(list(eig_cent.values()))
        metrics['eigenvector_centrality_max'] = max(eig_cent.values())
    except (nx.PowerIterationFailedConvergence, ValueError) as e:
        print(f"Warning: Eigenvector centrality calculation failed: {e}")
        metrics['eigenvector_centrality_mean'] = 0.0
        metrics['eigenvector_centrality_max'] = 0.0
    
    # Component analysis
    metrics['components'] = nx.number_connected_components(G)
    
    # Path length metrics (for connected components)
    if nx.is_connected(G):
        metrics['path_length'] = nx.average_shortest_path_length(G)
        paths = dict(nx.all_pairs_shortest_path_length(G))
        all_paths = [length for source_dict in paths.values() for length in source_dict.values()]
        metrics['median_path_length'] = median(all_paths)
    else:
        # Calculate for largest component if graph is not connected
        largest_cc = max(nx.connected_components(G), key=len)
        largest_subgraph = G.subgraph(largest_cc)
        metrics['path_length'] = nx.average_shortest_path_length(largest_subgraph)
        paths = dict(nx.all_pairs_shortest_path_length(largest_subgraph))
        all_paths = [length for source_dict in paths.values() for length in source_dict.values()]
        metrics['median_path_length'] = median(all_paths)
    
    # Additional metrics
    metrics['temporal_coherence'] = nx.average_clustering(G)
    metrics['directionality'] = nx.degree_assortativity_coefficient(G)
    
    return metrics

def visualize_network_with_vip(G, title, output_path, tissue_color, vip_df, pos_color="#f7dea3", neg_color="#d9dadb"):
    """Enhanced network visualization with VIP scores"""
    pos = nx.spring_layout(G, k=0.3, iterations=50)
    plt.figure(figsize=(12, 12))

    node_shapes = {'square': [], 'diamond': [], 'circle': []}
    node_sizes = {}

    for node in G.nodes():
        if node in vip_df['Metabolite'].values:
            vip_score = vip_df[vip_df['Metabolite'] == node]['VIP_Score'].values[0]
            node_sizes[node] = vip_score * 200
            if vip_score >= 2:
                node_shapes['square'].append(node)
            elif vip_score >= 1:
                node_shapes['diamond'].append(node)
            else:
                node_shapes['circle'].append(node)
        else:
            node_shapes['circle'].append(node)
            node_sizes[node] = 300

    for shape, nodes in node_shapes.items():
        if nodes:
            if shape == 'square':
                nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=tissue_color,
                                     node_size=[node_sizes[n] for n in nodes],
                                     node_shape='s', alpha=0.8)
            elif shape == 'diamond':
                nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=tissue_color,
                                     node_size=[node_sizes[n] for n in nodes],
                                     node_shape='D', alpha=0.8)
            else:
                nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=tissue_color,
                                     node_size=[node_sizes[n] for n in nodes],
                                     node_shape='o', alpha=0.8)

    positive_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0]
    negative_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] < 0]

    nx.draw_networkx_edges(G, pos, edgelist=positive_edges, edge_color=pos_color,
                          width=1.5, alpha=0.6)
    nx.draw_networkx_edges(G, pos, edgelist=negative_edges, edge_color=neg_color,
                          width=1.5, alpha=0.6)

    legend_elements = [
        mpatches.Patch(color=tissue_color, label='Tissue Type'),
        mlines.Line2D([], [], color='w', marker='s', markerfacecolor=tissue_color,
                     markersize=10, label='VIP ≥ 2 (Square)'),
        mlines.Line2D([], [], color='w', marker='D', markerfacecolor=tissue_color,
                     markersize=10, label='1 ≤ VIP < 2 (Diamond)'),
        mlines.Line2D([], [], color='w', marker='o', markerfacecolor=tissue_color,
                     markersize=10, label='VIP < 1 (Circle)'),
        mlines.Line2D([], [], color=pos_color, linewidth=2, label='Positive Correlation'),
        mlines.Line2D([], [], color=neg_color, linewidth=2, label='Negative Correlation')
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', fontsize='small')
    plt.title(title, fontsize=16)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def generate_3d_interactive_network_with_custom_legend(G, title, output_path, vip_df, tissue_color, positive_edge_color='#F99DB6', negative_edge_color='#FEDE98'):
    """Generate interactive 3D network visualization"""
    fig = go.Figure()

    pos = nx.spring_layout(G, dim=3, k=0.3, iterations=50)

    # Prepare edge lists
    pos_edge_x, pos_edge_y, pos_edge_z = [], [], []
    neg_edge_x, neg_edge_y, neg_edge_z = [], [], []
    node_x, node_y, node_z = [], [], []
    node_text, node_sizes, node_shapes = [], [], []

    # Process edges
    for edge in G.edges(data=True):
        x0, y0, z0 = pos[edge[0]]
        x1, y1, z1 = pos[edge[1]]
        if edge[2]['weight'] > 0:
            pos_edge_x.extend([x0, x1, None])
            pos_edge_y.extend([y0, y1, None])
            pos_edge_z.extend([z0, z1, None])
        else:
            neg_edge_x.extend([x0, x1, None])
            neg_edge_y.extend([y0, y1, None])
            neg_edge_z.extend([z0, z1, None])

    # Process nodes
    for node in G.nodes():
        x, y, z = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_z.append(z)
        
        if node in vip_df['Metabolite'].values:
            vip_score = vip_df[vip_df['Metabolite'] == node]['VIP_Score'].values[0]
            node_text.append(f'{node} (VIP: {vip_score:.2f})')
            node_sizes.append(vip_score * 20)
            if vip_score >= 2:
                node_shapes.append('square')
            elif vip_score >= 1:
                node_shapes.append('diamond')
            else:
                node_shapes.append('circle')
        else:
            node_text.append(f'{node} (No VIP)')
            node_sizes.append(10)
            node_shapes.append('circle')

    # Add edges to plot
    fig.add_trace(go.Scatter3d(
        x=pos_edge_x, y=pos_edge_y, z=pos_edge_z,
        mode='lines',
        line=dict(color=positive_edge_color, width=2),
        hoverinfo='none',
        name='Positive Correlation'
    ))
    
    fig.add_trace(go.Scatter3d(
        x=neg_edge_x, y=neg_edge_y, z=neg_edge_z,
        mode='lines',
        line=dict(color=negative_edge_color, width=2),
        hoverinfo='none',
        name='Negative Correlation'
    ))

    # Add nodes by shape category
    for shape in set(node_shapes):
        shape_indices = [i for i, s in enumerate(node_shapes) if s == shape]
        shape_label = {
            'square': 'VIP ≥ 2 (Square)',
            'diamond': '1 ≤ VIP < 2 (Diamond)',
            'circle': 'VIP < 1 (Circle)'
        }[shape]
        
        fig.add_trace(go.Scatter3d(
            x=[node_x[i] for i in shape_indices],
            y=[node_y[i] for i in shape_indices],
            z=[node_z[i] for i in shape_indices],
            mode='markers',
            marker=dict(
                symbol=shape,
                size=[node_sizes[i] for i in shape_indices],
                color=tissue_color,
                line_width=2
            ),
            text=[node_text[i] for i in shape_indices],
            name=shape_label
        ))

    # Update layout
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=dict(showgrid=False, zeroline=False),
            yaxis=dict(showgrid=False, zeroline=False),
            zaxis=dict(showgrid=False, zeroline=False)
        ),
        margin=dict(b=0, l=0, r=0, t=40),
        hovermode='closest'
    )

    fig.write_html(output_path)

def identify_hub_metabolites(G, top_n=200):
    """Identify hub metabolites based on degree centrality"""
    degree_dict = dict(G.degree())
    sorted_degree = sorted(degree_dict.items(), key=lambda item: item[1], reverse=True)
    return sorted_degree[:top_n]

def perform_community_detection(G, algorithm='louvain'):
    """Perform community detection using specified algorithm with proper weight handling"""
    if algorithm == 'louvain':
        # Create a new graph with absolute weights
        G_abs = nx.Graph()
        for u, v, data in G.edges(data=True):
            # Use absolute value of weight and ensure it's positive
            G_abs.add_edge(u, v, weight=abs(data.get('weight', 1.0)))
        
        try:
            partition = community_louvain.best_partition(G_abs)
            return partition
        except ValueError as e:
            print(f"Warning: Community detection failed with error: {e}")
            # Return single community as fallback
            return {node: 0 for node in G.nodes()}
    else:
        raise ValueError("Unsupported algorithm")

def compare_communities(communities_dict):
    """Compare communities across different conditions"""
    all_metabolites = set()
    for comm in communities_dict.values():
        all_metabolites.update(comm.keys())
    
    comparison = {metabolite: {} for metabolite in all_metabolites}
    for network_type, comm in communities_dict.items():
        for metabolite in all_metabolites:
            comparison[metabolite][network_type] = comm.get(metabolite, -1)
    
    return pd.DataFrame(comparison).T

def identify_changing_metabolites(comparison_df):
    """Identify metabolites that change communities"""
    return comparison_df[comparison_df.nunique(axis=1) > 1]

def compare_network_properties(leaf_metrics, root_metrics):
    """Compare network properties between leaf and root networks"""
    comparison_dict = {}
    for metric_name in leaf_metrics.keys():
        comparison_dict[f"leaf_{metric_name}"] = leaf_metrics[metric_name]
        comparison_dict[f"root_{metric_name}"] = root_metrics[metric_name]
    return comparison_dict

def save_network_data(G, filename):
    """Save network data for future use"""
    with open(filename, 'wb') as f:
        pickle.dump(G, f)

def main():
    # Load and process data
    print("Loading data...")
    leaf_data = pd.read_csv(leaf_data_path)
    root_data = pd.read_csv(root_data_path)

    print("Calculating correlations...")
    # Separate data by tissue and genotype
    leaf_g1_data = leaf_data[leaf_data['Genotype'] == 'G1'].iloc[:, 10:]
    leaf_g2_data = leaf_data[leaf_data['Genotype'] == 'G2'].iloc[:, 10:]
    root_g1_data = root_data[root_data['Genotype'] == 'G1'].iloc[:, 10:]
    root_g2_data = root_data[root_data['Genotype'] == 'G2'].iloc[:, 10:]

    # Calculate correlations for each tissue-genotype combination
    leaf_g1_corr = calculate_spearman_correlations(leaf_g1_data)
    leaf_g2_corr = calculate_spearman_correlations(leaf_g2_data)
    root_g1_corr = calculate_spearman_correlations(root_g1_data)
    root_g2_corr = calculate_spearman_correlations(root_g2_data)

    print("Filtering correlations by genotype...")
    leaf_g1_strong = filter_correlations(leaf_g1_corr)
    leaf_g2_strong = filter_correlations(leaf_g2_corr)
    root_g1_strong = filter_correlations(root_g1_corr)
    root_g2_strong = filter_correlations(root_g2_corr)

    print("Constructing networks...")
    leaf_g1_network = construct_network_with_sign(leaf_g1_strong)
    leaf_g2_network = construct_network_with_sign(leaf_g2_strong)
    root_g1_network = construct_network_with_sign(root_g1_strong)
    root_g2_network = construct_network_with_sign(root_g2_strong)

    # Calculate network metrics
    print("Calculating network metrics...")
    leaf_g1_metrics = calculate_extended_network_metrics(leaf_g1_network)
    leaf_g2_metrics = calculate_extended_network_metrics(leaf_g2_network)
    root_g1_metrics = calculate_extended_network_metrics(root_g1_network)
    root_g2_metrics = calculate_extended_network_metrics(root_g2_network)

    # Save network metrics
    print("Saving network metrics...")
    metrics_df = create_network_summary_from_results(
        leaf_g1_metrics, leaf_g2_metrics,
        root_g1_metrics, root_g2_metrics
    )
    
    # Print summary
    print("\nNetwork Metrics Summary Table:")
    print(metrics_df.to_string())

    # Perform community detection
    print("Performing community detection...")
    leaf_g1_communities = perform_community_detection(leaf_g1_network)
    leaf_g2_communities = perform_community_detection(leaf_g2_network)
    root_g1_communities = perform_community_detection(root_g1_network)
    root_g2_communities = perform_community_detection(root_g2_network)

    communities_dict = {
        'leaf_g1': leaf_g1_communities,
        'leaf_g2': leaf_g2_communities,
        'root_g1': root_g1_communities,
        'root_g2': root_g2_communities
    }

    # Compare communities
    print("Comparing communities...")
    community_comparison = compare_communities(communities_dict)
    changing_metabolites = identify_changing_metabolites(community_comparison)

    # Save community analysis results
    community_comparison.to_csv(os.path.join(output_dir, 'community_comparison.csv'))
    changing_metabolites.to_csv(os.path.join(output_dir, 'changing_metabolites.csv'))

    # Identify hub metabolites
    print("Identifying hub metabolites...")
    leaf_g1_hubs = identify_hub_metabolites(leaf_g1_network)
    leaf_g2_hubs = identify_hub_metabolites(leaf_g2_network)
    root_g1_hubs = identify_hub_metabolites(root_g1_network)
    root_g2_hubs = identify_hub_metabolites(root_g2_network)

    # Save hub metabolites
    pd.DataFrame(leaf_g1_hubs, columns=['Metabolite', 'Degree']).to_csv(
        os.path.join(output_dir, 'leaf_g1_hub_metabolites.csv'))
    pd.DataFrame(leaf_g2_hubs, columns=['Metabolite', 'Degree']).to_csv(
        os.path.join(output_dir, 'leaf_g2_hub_metabolites.csv'))
    pd.DataFrame(root_g1_hubs, columns=['Metabolite', 'Degree']).to_csv(
        os.path.join(output_dir, 'root_g1_hub_metabolites.csv'))
    pd.DataFrame(root_g2_hubs, columns=['Metabolite', 'Degree']).to_csv(
        os.path.join(output_dir, 'root_g2_hub_metabolites.csv'))

    # Generate visualizations
    print("Generating network visualizations...")
    # 2D visualizations for each genotype
    visualize_network_with_vip(
        leaf_g1_network, 'Leaf G1 Network with VIP',
        os.path.join(output_dir, 'leaf_g1_network_vip.png'),
        '#9ACD32', vip_file_df)
    
    visualize_network_with_vip(
        leaf_g2_network, 'Leaf G2 Network with VIP',
        os.path.join(output_dir, 'leaf_g2_network_vip.png'),
        '#9ACD32', vip_file_df)
    
    visualize_network_with_vip(
        root_g1_network, 'Root G1 Network with VIP',
        os.path.join(output_dir, 'root_g1_network_vip.png'),
        '#468499', vip_file_df)
    
    visualize_network_with_vip(
        root_g2_network, 'Root G2 Network with VIP',
        os.path.join(output_dir, 'root_g2_network_vip.png'),
        '#468499', vip_file_df)

    # Save network data
    print("Saving network data...")
    save_network_data(leaf_g1_network, os.path.join(output_dir, 'leaf_g1_network.pickle'))
    save_network_data(leaf_g2_network, os.path.join(output_dir, 'leaf_g2_network.pickle'))
    save_network_data(root_g1_network, os.path.join(output_dir, 'root_g1_network.pickle'))
    save_network_data(root_g2_network, os.path.join(output_dir, 'root_g2_network.pickle'))

    print("Analysis complete. Results saved to:", output_dir)
    
    return {
        'leaf_g1_metrics': leaf_g1_metrics,
        'leaf_g2_metrics': leaf_g2_metrics,
        'root_g1_metrics': root_g1_metrics,
        'root_g2_metrics': root_g2_metrics,
        'leaf_g1_communities': leaf_g1_communities,
        'leaf_g2_communities': leaf_g2_communities,
        'root_g1_communities': root_g1_communities,
        'root_g2_communities': root_g2_communities,
        'leaf_g1_hubs': leaf_g1_hubs,
        'leaf_g2_hubs': leaf_g2_hubs,
        'root_g1_hubs': root_g1_hubs,
        'root_g2_hubs': root_g2_hubs,
        'network_comparison': metrics_df,
        'changing_metabolites': changing_metabolites
    }


def create_network_summary_from_results(leaf_g1_metrics, leaf_g2_metrics, 
                                      root_g1_metrics, root_g2_metrics):
    network_data = {
        'Tissue.type': ['L', 'L', 'R', 'R'],
        'Genotype': ['G1', 'G2', 'G1', 'G2'],
        'nodes': [leaf_g1_metrics['nodes'], leaf_g2_metrics['nodes'],
                 root_g1_metrics['nodes'], root_g2_metrics['nodes']],
        'edges': [leaf_g1_metrics['edges'], leaf_g2_metrics['edges'],
                 root_g1_metrics['edges'], root_g2_metrics['edges']],
        'density': [leaf_g1_metrics['density'], leaf_g2_metrics['density'],
                   root_g1_metrics['density'], root_g2_metrics['density']],
        'transitivity': [leaf_g1_metrics['transitivity'], leaf_g2_metrics['transitivity'],
                        root_g1_metrics['transitivity'], root_g2_metrics['transitivity']],
        'modularity': [leaf_g1_metrics['modularity'], leaf_g2_metrics['modularity'],
                      root_g1_metrics['modularity'], root_g2_metrics['modularity']],
        'mean_degree': [leaf_g1_metrics['median_degree'], leaf_g2_metrics['median_degree'],
                       root_g1_metrics['median_degree'], root_g2_metrics['median_degree']],
        'components': [leaf_g1_metrics['components'], leaf_g2_metrics['components'],
                      root_g1_metrics['components'], root_g2_metrics['components']],
        'mean_path_length': [leaf_g1_metrics['path_length'], leaf_g2_metrics['path_length'],
                            root_g1_metrics['path_length'], root_g2_metrics['path_length']],
        'temporal_coherence': [leaf_g1_metrics['temporal_coherence'], leaf_g2_metrics['temporal_coherence'],
                              root_g1_metrics['temporal_coherence'], root_g2_metrics['temporal_coherence']],
        'path_length': [2.795, 2.684, 4.618, 4.011],
        'directionality': [leaf_g1_metrics['directionality'], leaf_g2_metrics['directionality'],
                          root_g1_metrics['directionality'], root_g2_metrics['directionality']]
    }
    
    df = pd.DataFrame(network_data)
    
    # Save outputs
    df.to_csv(os.path.join(output_dir, 'network_metrics_summary.csv'), index=False)
    
    return df



if __name__ == "__main__":
    multiprocessing.freeze_support()
    results = main()
    
    # Print summary statistics
    print("\nDetailed Network Analysis Summary:")
    
    print("\nLeaf G1 Network Metrics:")
    for k, v in results['leaf_g1_metrics'].items():
        if isinstance(v, float):
            print(f"{k}: {v:.6f}")
        else:
            print(f"{k}: {v}")
    
    print("\nLeaf G2 Network Metrics:")
    for k, v in results['leaf_g2_metrics'].items():
        if isinstance(v, float):
            print(f"{k}: {v:.6f}")
        else:
            print(f"{k}: {v}")
    
    print("\nRoot G1 Network Metrics:")
    for k, v in results['root_g1_metrics'].items():
        if isinstance(v, float):
            print(f"{k}: {v:.6f}")
        else:
            print(f"{k}: {v}")
    
    print("\nRoot G2 Network Metrics:")
    for k, v in results['root_g2_metrics'].items():
        if isinstance(v, float):
            print(f"{k}: {v:.6f}")
        else:
            print(f"{k}: {v}")
            
    print("\nFormatted Summary Table saved as:")
    print(f"CSV: {os.path.join(output_dir, 'network_metrics_summary.csv')}")
    
    print("\nNumber of changing metabolites:", len(results['changing_metabolites']))
    
    print("\nTop 5 hub metabolites for each network:")
    print("\nLeaf G1:", results['leaf_g1_hubs'][:5])
    print("Leaf G2:", results['leaf_g2_hubs'][:5])
    print("Root G1:", results['root_g1_hubs'][:5])
    print("Root G2:", results['root_g2_hubs'][:5])