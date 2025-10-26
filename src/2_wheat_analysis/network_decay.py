import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import mannwhitneyu, spearmanr, ks_2samp
import plotly.graph_objects as go
from statsmodels.nonparametric.smoothers_lowess import lowess
import os

# Set file paths
input_dir = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\spearman\network2_plot4E'
output_dir = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\spearman\network2_plot4E\decay'

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read the CSV files
leaf_g1 = pd.read_csv(os.path.join(input_dir, 'leaf_g1_hub_metabolites.csv'))
leaf_g2 = pd.read_csv(os.path.join(input_dir, 'leaf_g2_hub_metabolites.csv'))
root_g1 = pd.read_csv(os.path.join(input_dir, 'root_g1_hub_metabolites.csv'))
root_g2 = pd.read_csv(os.path.join(input_dir, 'root_g2_hub_metabolites.csv'))

# Add tissue and genotype labels
leaf_g1['group'] = 'Leaf_G1'
leaf_g2['group'] = 'Leaf_G2'
root_g1['group'] = 'Root_G1'
root_g2['group'] = 'Root_G2'

# Combine data
combined_df = pd.concat([leaf_g1, leaf_g2, root_g1, root_g2])

# 1. Kruskal-Wallis test
groups = [leaf_g1['Degree'], leaf_g2['Degree'], 
         root_g1['Degree'], root_g2['Degree']]
h_stat, p_value = stats.kruskal(*groups)

results = {
    'Kruskal-Wallis Test': {
        'H-statistic': h_stat,
        'p-value': p_value
    }
}

# 2. Mann-Whitney U tests for pairwise comparisons
pairs = [
    ('Leaf G1 vs G2', leaf_g1['Degree'], leaf_g2['Degree']),
    ('Root G1 vs G2', root_g1['Degree'], root_g2['Degree']),
    ('G1 Leaf vs Root', leaf_g1['Degree'], root_g1['Degree']),
    ('G2 Leaf vs Root', leaf_g2['Degree'], root_g2['Degree'])
]

mann_whitney_results = {}
for name, group1, group2 in pairs:
    stat, p_val = mannwhitneyu(group1, group2, alternative='two-sided')
    mann_whitney_results[name] = {'statistic': stat, 'p-value': p_val}
results['Mann-Whitney Tests'] = mann_whitney_results

# 3. Spearman correlations between degree and rank
spearman_results = {}
for name, df in [('Leaf_G1', leaf_g1), ('Leaf_G2', leaf_g2), 
                 ('Root_G1', root_g1), ('Root_G2', root_g2)]:
    df['rank'] = range(1, len(df) + 1)
    corr, p_val = spearmanr(df['rank'], df['Degree'])
    spearman_results[name] = {'correlation': corr, 'p-value': p_val}
results['Spearman Correlations'] = spearman_results

# 4. KS tests for distribution comparisons
ks_results = {}
for pair_name, group1, group2 in pairs:
    stat, p_val = ks_2samp(group1, group2)
    ks_results[pair_name] = {'statistic': stat, 'p-value': p_val}
results['Kolmogorov-Smirnov Tests'] = ks_results

# 5. LOWESS smoothing and visualization
fig = go.Figure()
colors = {'Leaf_G1': '#2ecc71', 'Leaf_G2': '#27ae60', 
          'Root_G1': '#33d6d3', 'Root_G2': '#1e597d'}

for name, df in [('Leaf_G1', leaf_g1), ('Leaf_G2', leaf_g2), 
                 ('Root_G1', root_g1), ('Root_G2', root_g2)]:
    # Raw data
    fig.add_trace(go.Scatter(
        x=df['rank'], y=df['Degree'],
        name=f'{name} (raw)',
        mode='markers',
        marker=dict(color=colors[name], size=5, opacity=0.5),
        showlegend=False
    ))
    
    # LOWESS smoothing
    smoothed = lowess(df['Degree'], df['rank'], frac=0.3)
    fig.add_trace(go.Scatter(
        x=smoothed[:, 0], y=smoothed[:, 1],
        name=name,
        line=dict(color=colors[name], width=2),
        mode='lines'
    ))

fig.update_layout(
    title='Hub Degree Decay Patterns with LOWESS Smoothing',
    xaxis_title='Hub Rank',
    yaxis_title='Degree',
    template='plotly_white'
)

# Save results
fig.write_html(os.path.join(output_dir, 'decay_patterns.html'))

# Save statistical results
with open(os.path.join(output_dir, 'statistical_analysis.txt'), 'w') as f:
    f.write("Network Decay Pattern Analysis Results\n")
    f.write("=====================================\n\n")
    
    f.write("1. Kruskal-Wallis Test\n")
    f.write(f"H-statistic: {results['Kruskal-Wallis Test']['H-statistic']:.4f}\n")
    f.write(f"p-value: {results['Kruskal-Wallis Test']['p-value']:.4e}\n\n")
    
    f.write("2. Mann-Whitney U Tests\n")
    for pair, result in results['Mann-Whitney Tests'].items():
        f.write(f"{pair}:\n")
        f.write(f"  Statistic: {result['statistic']:.4f}\n")
        f.write(f"  p-value: {result['p-value']:.4e}\n")
    f.write("\n")
    
    f.write("3. Spearman Correlations\n")
    for group, result in results['Spearman Correlations'].items():
        f.write(f"{group}:\n")
        f.write(f"  Correlation: {result['correlation']:.4f}\n")
        f.write(f"  p-value: {result['p-value']:.4e}\n")
    f.write("\n")
    
    f.write("4. Kolmogorov-Smirnov Tests\n")
    for pair, result in results['Kolmogorov-Smirnov Tests'].items():
        f.write(f"{pair}:\n")
        f.write(f"  Statistic: {result['statistic']:.4f}\n")
        f.write(f"  p-value: {result['p-value']:.4e}\n")

print("Analysis complete. Results saved to:", output_dir)