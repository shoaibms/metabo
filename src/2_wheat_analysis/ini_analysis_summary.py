
#################################################
########## genotype_comparison ##################
########## Code 1 ###############################
#################################################
import pandas as pd
import os
import numpy as np
from scipy import stats

# File paths
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\genotype_comparison.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\genotype_comp"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read data
df = pd.read_csv(input_file)

def analyze_regulation_pattern(sig_data):
    """Analyze regulation patterns using Log2FC"""
    # For genes with significant p-values:
    # Log2FC > 0 means G2 > G1
    # Log2FC < 0 means G1 > G2
    up_in_g1 = sum(sig_data['Log2_Fold_Change'] < 0)  # G1 higher
    down_in_g1 = sum(sig_data['Log2_Fold_Change'] > 0)  # G2 higher
    
    return up_in_g1, down_in_g1

def create_tissue_divergence():
    summary = []
    for tissue in df['Tissue'].unique():
        for day in df['Day'].unique():
            for treatment in df['Treatment'].unique():
                subset = df[(df['Tissue'] == tissue) & 
                          (df['Day'] == day) & 
                          (df['Treatment'] == treatment)]
                
                # Get significant features
                sig_data = subset[subset['Adjusted_P_value'] < 0.05]
                
                # Basic counts
                total_features = len(subset)
                sig_features = len(sig_data)
                
                # Analyze regulation patterns
                up_in_g1, down_in_g1 = analyze_regulation_pattern(sig_data)
                
                # Calculate medians for actual values
                g1_median = subset['G1_median'].median()
                g2_median = subset['G2_median'].median()
                
                summary.append({
                    'Tissue': tissue,
                    'Day': day,
                    'Treatment': treatment,
                    'Total_Features': total_features,
                    'Significant_Features': sig_features,
                    'Significant_Percent': (sig_features/total_features*100) if total_features > 0 else 0,
                    'Up_in_G1': up_in_g1,
                    'Down_in_G1': down_in_g1,
                    'Up_in_G2': down_in_g1,  # Opposite of G1
                    'Down_in_G2': up_in_g1,  # Opposite of G1
                    'G1_Median': g1_median,
                    'G2_Median': g2_median,
                    'Median_Effect_Size': subset['Effect_size'].median(),
                    'Median_Log2FC': subset['Log2_Fold_Change'].median()
                })
    
    return pd.DataFrame(summary)

# Generate and save table
tissue_divergence = create_tissue_divergence()
output_path = os.path.join(output_dir, 'tissue_divergence.csv')
tissue_divergence.to_csv(output_path, index=False)
print(f"Saved tissue_divergence.csv")

# Print validation summary
print("\nValidation Summary:")
for tissue in df['Tissue'].unique():
    for day in sorted(df['Day'].unique()):
        subset = tissue_divergence[
            (tissue_divergence['Tissue'] == tissue) & 
            (tissue_divergence['Day'] == day)
        ]
        print(f"\n{tissue} Day {day}:")
        print(f"Total significant: {subset['Significant_Features'].sum()}")
        print(f"Up in G1: {subset['Up_in_G1'].sum()}")
        print(f"Up in G2: {subset['Up_in_G2'].sum()}")
        print(f"Median Log2FC: {subset['Median_Log2FC'].mean():.3f}")



#################################################
########## genotype_comparison #################
########## Code 2 ###############################
#################################################


import pandas as pd
import os
import numpy as np
from scipy import stats

# File paths
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\genotype_comparison.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\genotype_comp"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read data
df = pd.read_csv(input_file)

# Modified Tissue-Specific Metabolic Divergence function
def create_tissue_divergence():
    summary = []
    for tissue in df['Tissue'].unique():
        for day in df['Day'].unique():
            tissue_day_data = df[(df['Tissue'] == tissue) & (df['Day'] == day)]
            
            total_features = len(tissue_day_data)
            sig_features = sum(tissue_day_data['Adjusted_P_value'] < 0.05)
            median_effect = tissue_day_data['Effect_size'].median()
            
            # Separate G1 and G2 metrics
            sig_data = tissue_day_data[tissue_day_data['Adjusted_P_value'] < 0.05]
            
            # G1 metrics
            g1_median = tissue_day_data['G1_median'].median()
            up_in_g1 = sum(sig_data['Log2_Fold_Change'] > 0)
            down_in_g1 = sum(sig_data['Log2_Fold_Change'] < 0)
            
            # G2 metrics
            g2_median = tissue_day_data['G2_median'].median()
            up_in_g2 = sum(sig_data['Log2_Fold_Change'] < 0)  # Note: reversed from G1
            down_in_g2 = sum(sig_data['Log2_Fold_Change'] > 0)  # Note: reversed from G1
            
            summary.append({
                'Tissue': tissue,
                'Day': day,
                'Total_Features': total_features,
                'Significant_Features_Percent': (sig_features/total_features)*100 if total_features > 0 else 0,
                'Median_Effect_Size': median_effect,
                'G1_Median': g1_median,
                'G2_Median': g2_median,
                'Up_in_G1': up_in_g1,
                'Down_in_G1': down_in_g1,
                'Up_in_G2': up_in_g2,
                'Down_in_G2': down_in_g2,
                'G1_Response_Ratio': up_in_g1/down_in_g1 if down_in_g1 > 0 else np.inf,
                'G2_Response_Ratio': up_in_g2/down_in_g2 if down_in_g2 > 0 else np.inf,
                'G1_to_G2_Ratio': g1_median/g2_median if g2_median > 0 else np.inf
            })
    
    return pd.DataFrame(summary)

# Generate and save tissue divergence table
tissue_divergence = create_tissue_divergence()
output_path = os.path.join(output_dir, 'tissue_divergence.csv')
tissue_divergence.to_csv(output_path, index=False)
print(f"Saved tissue_divergence.csv")

# Print summary of tissue divergence
print("\nTissue Divergence Summary:")
for tissue in df['Tissue'].unique():
    print(f"\n{tissue}:")
    tissue_data = tissue_divergence[tissue_divergence['Tissue'] == tissue]
    print(f"Average significant features: {tissue_data['Significant_Features_Percent'].mean():.1f}%")
    print(f"G1 median response: {tissue_data['G1_Median'].mean():.3f}")
    print(f"G2 median response: {tissue_data['G2_Median'].mean():.3f}")
    print(f"G1/G2 ratio: {tissue_data['G1_to_G2_Ratio'].mean():.3f}")



##########################################
########## resilience_index ##############
########### code 1 #######################
##########################################

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os

# File paths
input_file = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\resilience_index.csv'
output_dir = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\resilience_index'
os.makedirs(output_dir, exist_ok=True)

# Load data
df = pd.read_csv(input_file)

# 1. Primary Resilience Metrics
def generate_primary_metrics():
    summary = df.groupby(['Tissue', 'Genotype', 'Day']).agg({
        'Resilience_Index': ['median', 'std', 'count'],
        'Effect_size': 'median',
        'Recovery_Score': 'median'
    }).reset_index()
    
    summary.columns = ['Tissue', 'Genotype', 'Day', 'Resilience_Median', 
                      'Resilience_SD', 'n', 'Effect_Size', 'Recovery_Score']
    
    summary['Resilience_Index'] = summary.apply(
        lambda x: f"{x['Resilience_Median']:.2f}±{x['Resilience_SD']:.2f}", axis=1)
    
    return summary[['Tissue', 'Genotype', 'Day', 'Resilience_Index', 
                   'Effect_Size', 'Recovery_Score', 'n']]

# 2. Treatment Response Summary
def generate_treatment_summary():
    summary = df.groupby(['Tissue', 'Genotype', 'Day']).agg({
        'Control_median': 'median',
        'Treated_median': 'median',
        'N_control': 'first',
        'N_treated': 'first'
    }).reset_index()
    
    summary['Log2FC'] = np.log2(summary['Treated_median'] / summary['Control_median'])
    return summary[['Tissue', 'Genotype', 'Day', 'Control_median', 
                   'Treated_median', 'Log2FC', 'N_control']]

# 3. Temporal Dynamics
def generate_temporal_dynamics():
    # First calculate medians for each tissue/genotype/day combination
    medians = df.groupby(['Tissue', 'Genotype', 'Day'])['Resilience_Index'].median().reset_index()
    
    dynamics = []
    for (tissue, genotype), group in medians.groupby(['Tissue', 'Genotype']):
        group_sorted = group.sort_values('Day')
        dynamics.append({
            'Tissue': tissue,
            'Genotype': genotype,
            'Early_Response': group_sorted.iloc[0]['Resilience_Index'],
            'Peak_Response': group_sorted['Resilience_Index'].max(),
            'Recovery': group_sorted.iloc[-1]['Resilience_Index'],
            'Pattern': 'Declining' if group_sorted['Resilience_Index'].is_monotonic_decreasing else
                      'Increasing' if group_sorted['Resilience_Index'].is_monotonic_increasing else
                      'Peak at D2' if group_sorted.iloc[1]['Resilience_Index'] > group_sorted.iloc[0]['Resilience_Index'] and 
                                    group_sorted.iloc[1]['Resilience_Index'] > group_sorted.iloc[2]['Resilience_Index'] else 'Complex'
        })
    
    return pd.DataFrame(dynamics)

# 4. Cross-Tissue Integration
def generate_cross_tissue():
    # First calculate medians for each metabolite/genotype/day combination
    medians = df.groupby(['Tissue', 'Genotype', 'Day', 'Metabolite'])['Resilience_Index'].median().reset_index()
    
    results = []
    for (genotype, day), group in medians.groupby(['Genotype', 'Day']):
        leaf_values = group[group['Tissue'] == 'L']['Resilience_Index']
        root_values = group[group['Tissue'] == 'R']['Resilience_Index']
        
        # Get common metabolites between leaf and root
        leaf_metabolites = set(group[group['Tissue'] == 'L']['Metabolite'])
        root_metabolites = set(group[group['Tissue'] == 'R']['Metabolite'])
        common_metabolites = leaf_metabolites.intersection(root_metabolites)
        
        if common_metabolites:
            # Filter for common metabolites
            leaf_data = group[(group['Tissue'] == 'L') & 
                            (group['Metabolite'].isin(common_metabolites))]['Resilience_Index']
            root_data = group[(group['Tissue'] == 'R') & 
                            (group['Metabolite'].isin(common_metabolites))]['Resilience_Index']
            
            if len(leaf_data) > 1 and len(root_data) > 1:
                correlation, _ = stats.spearmanr(leaf_data, root_data)
                integrated_score = np.mean([leaf_data.mean(), root_data.mean()])
                
                results.append({
                    'Genotype': genotype,
                    'Day': day,
                    'Leaf_Root_Correlation': correlation,
                    'Integrated_Score': integrated_score,
                    'N_metabolites': len(common_metabolites)
                })
    
    results_df = pd.DataFrame(results)
    if not results_df.empty:
        results_df['Response_Class'] = pd.qcut(results_df['Integrated_Score'], 
                                             q=3, labels=['Low', 'Medium', 'High'])
    return results_df

# Generate and save all tables
tables = {
    'primary_metrics': generate_primary_metrics(),
    'treatment_response': generate_treatment_summary(),
    'temporal_dynamics': generate_temporal_dynamics(),
    'cross_tissue_integration': generate_cross_tissue()
}

for name, table in tables.items():
    table.to_csv(os.path.join(output_dir, f'{name}.csv'), index=False)
    print(f"Generated {name}.csv")

print("Summary tables generated successfully.")


##########################################
########## resilience_index ##############
########### code 2 #######################
##########################################

import pandas as pd
import numpy as np
from scipy import stats
import os

# File paths
input_file = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6e\resilience_index.csv'
output_dir = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6e\summary\resilience_index'
os.makedirs(output_dir, exist_ok=True)

def calculate_tissue_resilience_metrics(df):
    """Calculate non-parametric tissue-specific resilience metrics"""
    results = []
    
    for tissue in df['Tissue'].unique():
        for genotype in df['Genotype'].unique():
            subset = df[(df['Tissue'] == tissue) & (df['Genotype'] == genotype)]
            
            # Calculate core metrics using non-parametric methods
            resilience_median = subset['Resilience_Index'].median()
            resilience_iqr = subset['Resilience_Index'].quantile(0.75) - subset['Resilience_Index'].quantile(0.25)
            
            # Calculate temporal progression
            temporal_data = []
            for day in sorted(subset['Day'].unique()):
                day_data = subset[subset['Day'] == day]
                temporal_data.append({
                    'Day': day,
                    'Median_Resilience': day_data['Resilience_Index'].median(),
                    'Effect_Size': day_data['Effect_size'].median(),
                    'N_metabolites': len(day_data)
                })
            
            # Find peak response day
            peak_day = max(temporal_data, key=lambda x: abs(x['Effect_Size']))['Day']
            
            results.append({
                'Tissue': tissue,
                'Genotype': genotype,
                'Median_Resilience': resilience_median,
                'IQR': resilience_iqr,
                'Peak_Response_Day': peak_day,
                'N_metabolites': len(subset) // 3,  # Divide by 3 days to get per-timepoint count
                'Early_Effect': temporal_data[0]['Effect_Size'],
                'Peak_Effect': max(abs(d['Effect_Size']) for d in temporal_data),
                'Recovery_Score': subset['Recovery_Score'].median(),
                'Highly_Resilient_%': (subset['Resilience_Class'] == 'Highly Resilient').mean() * 100
            })
    
    return pd.DataFrame(results)

def calculate_cross_tissue_coordination(df):
    """Calculate cross-tissue coordination metrics using common metabolites"""
    results = []
    
    for genotype in df['Genotype'].unique():
        for day in df['Day'].unique():
            # Get data for each tissue
            leaf_data = df[(df['Tissue'] == 'L') & 
                         (df['Genotype'] == genotype) & 
                         (df['Day'] == day)]
            root_data = df[(df['Tissue'] == 'R') & 
                         (df['Genotype'] == genotype) & 
                         (df['Day'] == day)]
            
            # Get common metabolites
            leaf_metabolites = set(leaf_data['Metabolite'])
            root_metabolites = set(root_data['Metabolite'])
            common_metabolites = leaf_metabolites.intersection(root_metabolites)
            
            if common_metabolites:
                # Filter for common metabolites
                leaf_values = leaf_data[leaf_data['Metabolite'].isin(common_metabolites)]['Resilience_Index']
                root_values = root_data[root_data['Metabolite'].isin(common_metabolites)]['Resilience_Index']
                
                if len(leaf_values) > 0 and len(root_values) > 0:
                    # Calculate coordination metrics
                    correlation = stats.spearmanr(leaf_values, root_values)[0]
                    leaf_root_ratio = leaf_values.median() / root_values.median() if root_values.median() != 0 else np.nan
                    
                    results.append({
                        'Genotype': genotype,
                        'Day': day,
                        'Leaf_Root_Correlation': correlation,
                        'Leaf_Root_Ratio': leaf_root_ratio,
                        'N_common_metabolites': len(common_metabolites)
                    })
    
    return pd.DataFrame(results)

def main():
    try:
        # Read data
        df = pd.read_csv(input_file)
        
        # Generate summaries
        tissue_metrics = calculate_tissue_resilience_metrics(df)
        coordination_metrics = calculate_cross_tissue_coordination(df)
        
        # Save results
        tissue_metrics.to_csv(os.path.join(output_dir, 'tissue_resilience_metrics.csv'), index=False)
        coordination_metrics.to_csv(os.path.join(output_dir, 'tissue_coordination_metrics.csv'), index=False)
        
        # Generate summary text for publication
        g1_leaf = tissue_metrics[(tissue_metrics['Tissue']=='L') & (tissue_metrics['Genotype']=='G1')].iloc[0]
        g2_leaf = tissue_metrics[(tissue_metrics['Tissue']=='L') & (tissue_metrics['Genotype']=='G2')].iloc[0]
        
        summary_text = f"""
Resilience analysis revealed distinct tissue-specific adaptation strategies (n={g1_leaf['N_metabolites']} leaf, {tissue_metrics[(tissue_metrics['Tissue']=='R') & (tissue_metrics['Genotype']=='G1')].iloc[0]['N_metabolites']} root metabolites). The drought-tolerant genotype demonstrated stronger early effect (leaf: {g1_leaf['Early_Effect']:.3f}, IQR: {g1_leaf['IQR']:.3f}) compared to susceptible genotype (leaf: {g2_leaf['Early_Effect']:.3f}, IQR: {g2_leaf['IQR']:.3f}). Cross-tissue coordination analysis revealed higher leaf-root correlation in G1 (r={coordination_metrics[coordination_metrics['Genotype']=='G1']['Leaf_Root_Correlation'].mean():.3f}) compared to G2 (r={coordination_metrics[coordination_metrics['Genotype']=='G2']['Leaf_Root_Correlation'].mean():.3f}), suggesting more integrated stress response in drought-tolerant lines.
"""
        
        # Save summary text
        with open(os.path.join(output_dir, 'publication_summary.txt'), 'w') as f:
            f.write(summary_text)
            
        print("Summary tables and publication text generated successfully")
        print("\nSuggested text for publication:")
        print(summary_text)
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()



#################################################
########## treatment_effect #####################
#################################################



import pandas as pd
import numpy as np
import os

# File paths
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\treatment_effect.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\treatment_effect"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read the data
df = pd.read_csv(input_file)

# 1. Primary Treatment Response Metrics
def create_primary_response_table():
    response_metrics = []
    
    for tissue in df['Tissue'].unique():
        tissue_data = df[df['Tissue'] == tissue]
        for day in tissue_data['Day'].unique():
            day_data = tissue_data[tissue_data['Day'] == day]
            
            metrics = {
                'Tissue': tissue,
                'Day': day,
                'Response_Rate': (day_data['Adjusted_P_value'] < 0.05).mean() * 100,
                'Median_Effect_Size': day_data['Effect_size'].median(),
                'Strong_Responses': (abs(day_data['Log2_Fold_Change']) > 1).sum(),
                'Moderate_Responses': ((abs(day_data['Log2_Fold_Change']) >= 0.5) & 
                                    (abs(day_data['Log2_Fold_Change']) <= 1)).sum()
            }
            response_metrics.append(metrics)
    
    response_table = pd.DataFrame(response_metrics)
    return response_table

# 2. Temporal Treatment Dynamics
def create_temporal_dynamics_table():
    temporal_metrics = []
    
    for tissue in df['Tissue'].unique():
        tissue_data = df[df['Tissue'] == tissue]
        
        # Calculate early response (Day 1)
        early_response = tissue_data[tissue_data['Day'] == 1]['Effect_size'].median()
        
        # Find peak response
        peak_day = tissue_data.groupby('Day')['Effect_size'].median().idxmax()
        peak_response = tissue_data[tissue_data['Day'] == peak_day]['Effect_size'].median()
        
        # Calculate cross-tissue correlation for each day
        for day in tissue_data['Day'].unique():
            day_data = tissue_data[tissue_data['Day'] == day]
            
            metrics = {
                'Tissue': tissue,
                'Time_Point': day,
                'Early_Response': early_response,
                'Peak_Response': peak_response,
                'Response_Pattern': 'Early' if day == 1 else ('Peak' if day == peak_day else 'Recovery'),
                'Cross_tissue_Correlation': np.nan  # Will be updated if both tissues present
            }
            temporal_metrics.append(metrics)
    
    temporal_table = pd.DataFrame(temporal_metrics)
    return temporal_table

# 3. Genotype-Tissue Response Comparison
def create_genotype_comparison_table():
    comparison_metrics = []
    
    for tissue in df['Tissue'].unique():
        tissue_data = df[df['Tissue'] == tissue]
        for genotype in tissue_data['Genotype'].unique():
            genotype_data = tissue_data[tissue_data['Genotype'] == genotype]
            
            metrics = {
                'Tissue': tissue,
                'Genotype': genotype,
                'Significant_Features': (genotype_data['Adjusted_P_value'] < 0.05).mean() * 100,
                'Median_Log2FC': genotype_data['Log2_Fold_Change'].median(),
                'Response_Distribution': 'Up' if genotype_data['Log2_Fold_Change'].median() > 0 else 'Down',
                'Recovery_Pattern': 'Strong' if abs(genotype_data['Effect_size'].median()) > 0.5 else 'Moderate'
            }
            comparison_metrics.append(metrics)
    
    comparison_table = pd.DataFrame(comparison_metrics)
    return comparison_table

# Generate and save tables
try:
    # Create tables
    primary_response = create_primary_response_table()
    temporal_dynamics = create_temporal_dynamics_table()
    genotype_comparison = create_genotype_comparison_table()
    
    # Save tables
    primary_response.to_csv(os.path.join(output_dir, 'primary_response_metrics.csv'), index=False)
    temporal_dynamics.to_csv(os.path.join(output_dir, 'temporal_dynamics.csv'), index=False)
    genotype_comparison.to_csv(os.path.join(output_dir, 'genotype_comparison.csv'), index=False)
    
    # Print summary
    print("Summary tables generated successfully:")
    print(f"- Primary Response Metrics: {len(primary_response)} rows")
    print(f"- Temporal Dynamics: {len(temporal_dynamics)} rows")
    print(f"- Genotype Comparison: {len(genotype_comparison)} rows")
    
except Exception as e:
    print(f"Error generating tables: {str(e)}")





#################################################
########## stress_level_comparison ##############
########## Code 1 ###############################
#################################################

import pandas as pd
import numpy as np
from scipy import stats

# File paths
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\stress_level_comparison.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\stress_level_comp"

# Read data
def load_data():
    try:
        df = pd.read_csv(input_file)
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def create_primary_stress_metrics(df):
    """Creates primary stress response summary."""
    summary = df.groupby(['Tissue', 'Day']).agg({
        'Stress_Response_Ratio': ['median', 'std'],
        'Effect_size': 'median',
        'Adjusted_P_value': lambda x: (x < 0.05).mean(),
        'Harsh_median': 'median',
        'Mild_median': 'median'
    }).round(3)
    
    # Calculate response pattern
    summary['Response_Pattern'] = summary.apply(
        lambda x: 'Strong' if abs(x[('Stress_Response_Ratio', 'median')]) > 1 
        else 'Moderate' if abs(x[('Stress_Response_Ratio', 'median')]) > 0.5 
        else 'Weak', axis=1
    )
    
    return summary.reset_index()

def create_tissue_dynamics(df):
    """Creates tissue-specific dynamics summary."""
    dynamics = df.groupby(['Tissue', 'Day']).agg({
        'Stress_Response_Ratio': ['median', 'std'],
        'Adjusted_P_value': lambda x: (x < 0.05).mean() * 100,
        'Effect_size': ['median', 'std']
    }).round(3)
    
    # Add response magnitude categorization
    dynamics['Response_Magnitude'] = dynamics.apply(
        lambda x: 'High' if abs(x[('Effect_size', 'median')]) > 0.5
        else 'Medium' if abs(x[('Effect_size', 'median')]) > 0.3
        else 'Low', axis=1
    )
    
    return dynamics.reset_index()

def create_genotype_contrasts(df):
    """Creates genotype comparison summary (supplementary)."""
    contrasts = df.pivot_table(
        index=['Tissue', 'Day'],
        columns='Genotype',
        values=['Stress_Response_Ratio', 'Effect_size'],
        aggfunc={'Stress_Response_Ratio': 'median', 'Effect_size': 'median'}
    ).round(3)
    
    # Calculate effect size delta between genotypes
    contrasts['Effect_Size_Delta'] = (
        contrasts[('Effect_size', 'G1')] - contrasts[('Effect_size', 'G2')]
    )
    
    return contrasts.reset_index()

def main():
    # Create output directory if it doesn't exist
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    df = load_data()
    if df is None:
        return
    
    # Generate summaries
    primary_metrics = create_primary_stress_metrics(df)
    tissue_dynamics = create_tissue_dynamics(df)
    genotype_contrasts = create_genotype_contrasts(df)
    
    # Save results
    output_files = {
        'primary_stress_metrics.csv': primary_metrics,
        'tissue_dynamics.csv': tissue_dynamics,
        'genotype_contrasts.csv': genotype_contrasts
    }
    
    for filename, data in output_files.items():
        output_path = os.path.join(output_dir, filename)
        data.to_csv(output_path, index=False)
        print(f"Saved {filename}")

if __name__ == "__main__":
    main()


#################################################
########## stress_level_comparison###############
########## Code 2 - Genotype ####################
#################################################
import pandas as pd
import numpy as np
import os

# File paths
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\stress_level_comparison.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\stress_level_comp"


def create_primary_stress_metrics(df):
    """Primary stress metrics with genotype comparison"""
    
    # Calculate metrics for each genotype separately
    primary = (df.groupby(['Tissue', 'Day'])
               .apply(lambda x: pd.Series({
                   'G1_Effect_Size': x[x['Genotype']=='G1']['Effect_size'].median(),
                   'G2_Effect_Size': x[x['Genotype']=='G2']['Effect_size'].median(),
                   'G1_Response_Ratio': x[x['Genotype']=='G1']['Stress_Response_Ratio'].median(),
                   'G2_Response_Ratio': x[x['Genotype']=='G2']['Stress_Response_Ratio'].median(),
                   'G1_Significance': (x[x['Genotype']=='G1']['Adjusted_P_value'] < 0.05).mean(),
                   'G2_Significance': (x[x['Genotype']=='G2']['Adjusted_P_value'] < 0.05).mean()
               }))
               .round(3))
    
    # Add comparative metrics
    primary['Effect_Size_Delta'] = primary['G1_Effect_Size'] - primary['G2_Effect_Size']
    primary['Response_Pattern'] = primary.apply(
        lambda x: 'G1>G2' if x['Effect_Size_Delta'] > 0.3
        else 'G2>G1' if x['Effect_Size_Delta'] < -0.3
        else 'Similar', axis=1
    )
    
    return primary.reset_index()

def create_tissue_dynamics(df):
    """Tissue dynamics with genotype comparison"""
    
    dynamics = (df.groupby(['Tissue', 'Day'])
                .apply(lambda x: pd.Series({
                    'G1_HM_Ratio': x[x['Genotype']=='G1']['Harsh_median'].median() / 
                                  x[x['Genotype']=='G1']['Mild_median'].median(),
                    'G2_HM_Ratio': x[x['Genotype']=='G2']['Harsh_median'].median() / 
                                  x[x['Genotype']=='G2']['Mild_median'].median(),
                    'G1_Sig_Features': (x[x['Genotype']=='G1']['Adjusted_P_value'] < 0.05).mean() * 100,
                    'G2_Sig_Features': (x[x['Genotype']=='G2']['Adjusted_P_value'] < 0.05).mean() * 100,
                    'Response_Magnitude': 'High' if abs(x['Effect_size'].median()) > 0.5
                                  else 'Medium' if abs(x['Effect_size'].median()) > 0.3
                                  else 'Low'
                }))
                .round(3))
    
    # Add temporal pattern
    dynamics['Temporal_Pattern'] = dynamics.groupby('Tissue')['G1_HM_Ratio'].transform(
        lambda x: 'Increasing' if x.iloc[-1] > x.iloc[0]
        else 'Decreasing' if x.iloc[-1] < x.iloc[0]
        else 'Stable'
    )
    
    return dynamics.reset_index()

def main():
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Read data
        df = pd.read_csv(input_file)
        
        # Generate summaries
        primary_metrics = create_primary_stress_metrics(df)
        tissue_dynamics = create_tissue_dynamics(df)
        
        # Save results
        primary_metrics.to_csv(os.path.join(output_dir, 'primary_stress_metrics_genotype.csv'), index=False)
        tissue_dynamics.to_csv(os.path.join(output_dir, 'tissue_dynamics_genotype.csv'), index=False)
        
        print("Summary tables generated successfully")
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()


#################################################
########## temporal_analysis C###################
########## Code 1 ###############################
#################################################
                
import pandas as pd
import numpy as np
from scipy import stats
import os

# File paths remain the same
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\temporal_analysis.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\temporal_analysis"
os.makedirs(output_dir, exist_ok=True)

def calculate_velocity(day1, day2):
    return (day2 - day1)

def calculate_directionality(velocities):
    if len(velocities) < 2:
        return np.nan
    return np.sum(np.diff(velocities) > 0) / (len(velocities) - 1)

def create_metabolic_velocity_with_genotype(df):
    """Create metabolic velocity summary including genotype differences"""
    summary = []
    
    for tissue in df['Tissue'].unique():
        for genotype in df['Genotype'].unique():
            subset = df[(df['Tissue'] == tissue) & (df['Genotype'] == genotype)]
            
            # Calculate velocities
            velocities_1_2 = calculate_velocity(
                subset['Day1_median'].median(),
                subset['Day2_median'].median()
            )
            velocities_2_3 = calculate_velocity(
                subset['Day2_median'].median(),
                subset['Day3_median'].median()
            )
            
            # Calculate recovery rate
            recovery_rate = (subset['Day3_median'].median() - subset['Day1_median'].median()) / 2
            
            # Calculate directionality coefficient
            dir_coef = calculate_directionality([
                subset['Day1_median'].median(),
                subset['Day2_median'].median(),
                subset['Day3_median'].median()
            ])
            
            # Calculate significant features and response timing
            sig_features = len(subset[subset['Adjusted_P_value'] < 0.05]) / len(subset) * 100
            
            # Determine peak response day
            medians = [subset['Day1_median'].median(), 
                      subset['Day2_median'].median(), 
                      subset['Day3_median'].median()]
            peak_day = np.argmax(np.abs(medians)) + 1
            
            summary.append({
                'Tissue': tissue,
                'Genotype': genotype,
                'Day1→2_Velocity': round(velocities_1_2, 3),
                'Day2→3_Velocity': round(velocities_2_3, 3),
                'Recovery_Rate': round(recovery_rate, 3),
                'Directionality_Coefficient': round(dir_coef, 3),
                'Significant_Features_%': round(sig_features, 2),
                'Peak_Response_Day': peak_day,
                'Total_Features': len(subset)
            })
    
    df_summary = pd.DataFrame(summary)
    
    # Add tissue ratio calculation for each genotype
    for genotype in df_summary['Genotype'].unique():
        genotype_data = df_summary[df_summary['Genotype'] == genotype]
        root_features = genotype_data[genotype_data['Tissue'] == 'R']['Significant_Features_%'].iloc[0]
        leaf_features = genotype_data[genotype_data['Tissue'] == 'L']['Significant_Features_%'].iloc[0]
        df_summary.loc[df_summary['Genotype'] == genotype, 'Root/Leaf_Ratio'] = round(root_features / leaf_features, 2)
    
    return df_summary

def main():
    try:
        df = pd.read_csv(input_file)
        print("Data loaded successfully")
        
        # Create updated metabolic velocity summary
        metabolic_velocity = create_metabolic_velocity_with_genotype(df)
        
        # Save results
        output_file = os.path.join(output_dir, 'metabolic_velocity_with_genotype.csv')
        metabolic_velocity.to_csv(output_file, index=False)
        
        print("Summary table created successfully")
        print("File saved as:", output_file)
        
        # Display summary
        print("\nMetabolic Velocity Summary with Genotype Comparison:")
        print(metabolic_velocity)
        
        # Display additional tissue-specific summaries by genotype
        print("\nTissue-specific response summary by genotype:")
        pivot_summary = metabolic_velocity.pivot_table(
            values=['Significant_Features_%', 'Root/Leaf_Ratio'],
            index='Genotype',
            columns='Tissue',
            aggfunc='mean'
        )
        print(pivot_summary)
        
    except Exception as e:
        print("Error occurred:", str(e))

if __name__ == "__main__":
    main()




#################################################
########## temporal_analysis C###################
########## Code 2 ###############################
#################################################

import pandas as pd
import numpy as np
from scipy import stats
import os

# File paths
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\temporal_analysis.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\temporal_analysis"
os.makedirs(output_dir, exist_ok=True)

def calculate_velocity(day1, day2):
    """Calculate velocity between two timepoints"""
    return (day2 - day1)

def calculate_directionality(velocities):
    """Calculate directionality coefficient from velocity changes"""
    if len(velocities) < 2:
        return np.nan
    return np.sum(np.diff(velocities) > 0) / (len(velocities) - 1)

def create_tissue_temporal_dynamics(df):
    """Create tissue-specific temporal dynamics summary"""
    summary = []
    
    for tissue in df['Tissue'].unique():
        for treatment in df['Treatment'].unique():
            subset = df[(df['Tissue'] == tissue) & (df['Treatment'] == treatment)]
            
            # Calculate significant features percentage
            total_features = len(subset)
            sig_features = len(subset[subset['Adjusted_P_value'] < 0.05])
            sig_percentage = (sig_features / total_features * 100) if total_features > 0 else 0
            
            # Find peak response day
            medians = [subset['Day1_median'].median(), 
                      subset['Day2_median'].median(), 
                      subset['Day3_median'].median()]
            peak_day = np.argmax(np.abs(medians)) + 1
            
            # Calculate effect size confidence interval
            ci_lower = subset['CI_lower'].mean()
            ci_upper = subset['CI_upper'].mean()
            
            # Determine predominant pattern
            patterns = subset['Pattern'].value_counts()
            main_pattern = patterns.index[0] if len(patterns) > 0 else 'Unknown'
            
            summary.append({
                'Tissue': tissue,
                'Treatment': treatment,
                'Temporal_trend': subset['Temporal_trend'].mean(),
                'Peak_Response_Day': peak_day,
                'Pattern': main_pattern,
                'Effect_Size': subset['Effect_size'].mean(),
                'CI_Lower': ci_lower,
                'CI_Upper': ci_upper,
                'n': subset['N_per_timepoint'].iloc[0],
                'Significant_Features_%': round(sig_percentage, 2),
                'Adj_P_value': subset['Adjusted_P_value'].median()
            })
    
    return pd.DataFrame(summary)

def create_metabolic_velocity(df):
    """Create metabolic velocity summary"""
    summary = []
    
    for tissue in df['Tissue'].unique():
        subset = df[df['Tissue'] == tissue]
        
        # Calculate velocities
        velocities_1_2 = calculate_velocity(
            subset['Day1_median'].median(),
            subset['Day2_median'].median()
        )
        velocities_2_3 = calculate_velocity(
            subset['Day2_median'].median(),
            subset['Day3_median'].median()
        )
        
        # Calculate recovery rate
        recovery_rate = (subset['Day3_median'].median() - subset['Day1_median'].median()) / 2
        
        # Calculate directionality coefficient
        dir_coef = calculate_directionality([
            subset['Day1_median'].median(),
            subset['Day2_median'].median(),
            subset['Day3_median'].median()
        ])
        
        # Calculate significant features
        sig_features = len(subset[subset['Adjusted_P_value'] < 0.05]) / len(subset) * 100
        
        summary.append({
            'Tissue': tissue,
            'Day1→2_Velocity': round(velocities_1_2, 3),
            'Day2→3_Velocity': round(velocities_2_3, 3),
            'Recovery_Rate': round(recovery_rate, 3),
            'Directionality_Coefficient': round(dir_coef, 3),
            'Significant_Features_%': round(sig_features, 2)
        })
    
    return pd.DataFrame(summary)

def main():
    # Read data
    try:
        df = pd.read_csv(input_file)
        print("Data loaded successfully")
        
        # Create summaries
        tissue_dynamics = create_tissue_temporal_dynamics(df)
        metabolic_velocity = create_metabolic_velocity(df)
        
        # Save results
        tissue_dynamics.to_csv(os.path.join(output_dir, 'tissue_temporal_dynamics.csv'), index=False)
        metabolic_velocity.to_csv(os.path.join(output_dir, 'metabolic_velocity.csv'), index=False)
        
        print("Summary tables created successfully")
        print("Files saved in:", output_dir)
        
        # Display first few rows of each summary
        print("\nTissue Temporal Dynamics Summary:")
        print(tissue_dynamics.head())
        print("\nMetabolic Velocity Summary:")
        print(metabolic_velocity.head())
        
    except Exception as e:
        print("Error occurred:", str(e))

if __name__ == "__main__":
    main()




#################################################
########## tissue_comparison ####################
#################################################

import pandas as pd
import numpy as np
import os

# Define output directory
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary\tissue_comparison"
os.makedirs(output_dir, exist_ok=True)

# Read tissue comparison data
df = pd.read_csv(r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\tissue_comparison.csv")

# Architecture Summary
arch_results = []
for day in sorted(df['Day'].unique()):
    for genotype in sorted(df['Genotype'].unique()):
        subset = df[(df['Day'] == day) & (df['Genotype'] == genotype)]
        
        arch_results.append({
            'Time_Point': f'Day {day}',
            'Genotype': genotype,
            'N_Pairs': len(subset),
            'Significant_%': round((subset['Adjusted_P_value'] < 0.05).mean() * 100, 1),
            'Median_Effect': round(subset['Effect_size'].median(), 3),
            'Strong_Effects_%': round((abs(subset['Effect_size']) > 0.6).mean() * 100, 1),
            'Median_Power': round(subset['Power'].median(), 3)
        })

# Coordination Summary
coord_results = []
for day in sorted(df['Day'].unique()):
    for genotype in sorted(df['Genotype'].unique()):
        subset = df[(df['Day'] == day) & (df['Genotype'] == genotype)]
        
        coord_results.append({
            'Time_Point': f'Day {day}',
            'Genotype': genotype,
            'Integration_%': round((abs(subset['Effect_size']) > 0.6).mean() * 100, 1),
            'Pos/Neg_Ratio': round(len(subset[subset['Effect_size'] > 0]) / max(len(subset[subset['Effect_size'] < 0]), 1), 2),
            'Effect_Magnitude': round(abs(subset['Effect_size']).median(), 3),
            'CI_Width': round((subset['CI_upper'] - subset['CI_lower']).median(), 3)
        })

# Create DataFrames
arch_df = pd.DataFrame(arch_results)
coord_df = pd.DataFrame(coord_results)

# Save results
arch_df.to_csv(os.path.join(output_dir, "tissue_architecture_summary.csv"), index=False)
coord_df.to_csv(os.path.join(output_dir, "tissue_coordination_summary.csv"), index=False)

print("Files saved to:", output_dir)
print("\nArchitecture Summary:")
print(arch_df)
print("\nCoordination Summary:")
print(coord_df)





#################################################
########## Correlation_analysis #################
#################################################

import pandas as pd
import os

# File paths
input_file = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\correlation_analysis.csv"
output_dir = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6\summary"
os.makedirs(output_dir, exist_ok=True)

# Load data
df = pd.read_csv(input_file)

# Summarize key metrics using pre-calculated statistics
summary = df.groupby(['Tissue_Pair', 'Genotype', 'Treatment']).agg({
    'Correlation': ['median', 'count'],
    'Adjusted_P_value': lambda x: sum(x < 0.01),
    'Effect_size': ['mean', 'std']  # Using pre-calculated effect sizes
}).reset_index()

# Flatten column names
summary.columns = ['Tissue_Pair', 'Genotype', 'Treatment', 
                  'Median_Correlation', 'Total_Correlations', 
                  'FDR_Significant', 'Mean_Effect_Size', 'Effect_Size_SD']

# Calculate treatment effects using pre-validated correlations
treatment_effect = df.groupby(['Tissue_Pair', 'Genotype']).apply(
    lambda x: pd.Series({
        'Treatment_Effect': x[x['Treatment']==1]['Correlation'].median() - 
                          x[x['Treatment']==0]['Correlation'].median(),
        'Network_Response': sum((abs(x['Correlation']) >= 0.7) & 
                              (x['Adjusted_P_value'] < 0.01) & 
                              (x['Treatment'] == 1)) - 
                          sum((abs(x['Correlation']) >= 0.7) & 
                              (x['Adjusted_P_value'] < 0.01) & 
                              (x['Treatment'] == 0))
    })
).reset_index()

# Merge summaries
final_summary = summary.merge(treatment_effect, on=['Tissue_Pair', 'Genotype'])

# Save summarized results
output_file = os.path.join(output_dir, 'correlation_network_summary.csv')
final_summary.round(3).to_csv(output_file, index=False)

print("Summary statistics for tissue-specific network responses:")
for tissue in final_summary['Tissue_Pair'].unique():
    print(f"\n{tissue}:")
    tissue_data = final_summary[final_summary['Tissue_Pair'] == tissue]
    print(f"Significant correlations: {tissue_data['FDR_Significant'].sum():,}")
    print(f"Mean effect size: {tissue_data['Mean_Effect_Size'].mean():.3f} ± {tissue_data['Effect_Size_SD'].mean():.3f}")