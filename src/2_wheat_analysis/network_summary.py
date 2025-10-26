

####################################
###### hub_metabolites summary #####
####################################
import pandas as pd
import os

def analyze_hub_files(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # Load files
    files = {
        'leaf_g1': pd.read_csv(f'{input_dir}/leaf_g1_hub_metabolites.csv'),
        'leaf_g2': pd.read_csv(f'{input_dir}/leaf_g2_hub_metabolites.csv'),
        'root_g1': pd.read_csv(f'{input_dir}/root_g1_hub_metabolites.csv'), 
        'root_g2': pd.read_csv(f'{input_dir}/root_g2_hub_metabolites.csv')
    }
    
    # Calculate stats
    summary = pd.DataFrame({
        'tissue_genotype': ['Leaf G1', 'Leaf G2', 'Root G1', 'Root G2'],
        'max_degree': [files[f].Degree.max() for f in files],
        'min_degree': [files[f].Degree.min() for f in files],
        'mean_degree': [files[f].Degree.mean() for f in files],
        'top20_max': [files[f].Degree.iloc[:20].max() for f in files],
        'top20_min': [files[f].Degree.iloc[:20].min() for f in files],
        'top20_mean': [files[f].Degree.iloc[:20].mean() for f in files]
    })
    
    # Save summary
    summary.to_csv(f'{output_dir}/hub_summary.csv', index=False)
    
    # Count metabolite types
    type_summary = pd.DataFrame()
    for name, df in files.items():
        n_cluster = df.Metabolite.str.startswith('N_Cluster').sum()
        p_cluster = df.Metabolite.str.startswith('P_Cluster').sum()
        type_summary.loc[name, 'N_Cluster'] = n_cluster
        type_summary.loc[name, 'P_Cluster'] = p_cluster
    
    type_summary.to_csv(f'{output_dir}/metabolite_types.csv')

analyze_hub_files(
    'C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E',
    'C:/Users/ms/Desktop/data_chem_3_10/output/results/spearman/network2_plot4E/summary'
)



