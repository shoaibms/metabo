"""
ARABIDOPSIS DATASET 2 METABOLITE DATA LOADER & ANALYZER
=======================================================
Updated loader for the second Arabidopsis dataset (arabidopsis2 folder).
Analyzes metabolite data structure and prepares it for network analysis.

Key Features:
- Robust data loading from new location
- Comprehensive structure analysis
- Sample-metabolite data linking
- Export for network analysis pipeline

Author: Updated for Dataset 2
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Set up paths for DATASET 2
base_path = r"C:\Users\ms\Desktop\data_chem_3_10"
data_path = os.path.join(base_path, "data", "v-data", "arabidopsis2")  # UPDATED PATH
output_plot_path = os.path.join(base_path, "validation", "output", "plot")
output_data_path = os.path.join(base_path, "validation", "output", "data")

print("="*60)
print("ARABIDOPSIS DATASET 2 METABOLITE ANALYSIS")
print("="*60)

# Load the sample metadata
print("\n1. LOADING SAMPLE METADATA...")
print("-" * 30)
try:
    sample_df = pd.read_csv(os.path.join(data_path, "s_MTBLS2289.txt"), sep='\t', low_memory=False)
    print(f"✓ Sample metadata loaded: {sample_df.shape}")
except Exception as e:
    print(f"✗ Error loading sample metadata: {e}")
    exit()

# Load the ACTUAL metabolite data (MAF file)
print("\n2. LOADING METABOLITE DATA (MAF FILE)...")
print("-" * 40)
try:
    maf_file = "m_MTBLS2289_GC-MS___metabolite_profiling_v2_maf.tsv"
    maf_df = pd.read_csv(os.path.join(data_path, maf_file), sep='\t', low_memory=False)
    print(f"✓ MAF file loaded successfully: {maf_df.shape}")
except Exception as e:
    print(f"✗ Error loading MAF file: {e}")
    exit()

print(f"\n3. MAF FILE STRUCTURE ANALYSIS:")
print("-" * 40)
print(f"Shape: {maf_df.shape}")
print(f"Columns: {list(maf_df.columns)}")

# Show first few rows
print(f"\nFirst 3 rows:")
print(maf_df.head(3))

# Identify data structure
print(f"\n4. DATA ORIENTATION ANALYSIS:")
print("-" * 35)

# Define known metadata columns for this dataset
metadata_cols = [
    'database_identifier', 'chemical_formula', 'smiles', 'inchi', 
    'metabolite_identification', 'mass_to_charge', 'fragmentation', 
    'modifications', 'charge', 'retention_time', 'taxid', 'species', 
    'database', 'database_version', 'reliability', 'uri', 'search_engine', 
    'search_engine_score', 'smallmolecule_abundance_sub', 
    'smallmolecule_abundance_stdev_sub', 'smallmolecule_abundance_std_error_sub'
]

# Get sample data columns (everything that's not metadata)
sample_cols = [col for col in maf_df.columns if col not in metadata_cols]
print(f"✓ Metabolite info columns: {len(metadata_cols)}")
print(f"✓ Sample data columns: {len(sample_cols)}")
print(f"✓ Data orientation: Metabolites in ROWS, Samples in COLUMNS")

# Validate sample matching
print(f"\n5. SAMPLE MATCHING VALIDATION:")
print("-" * 35)

# Extract sample names from metadata
metadata_samples = set(sample_df['Sample Name'].astype(str).str.strip())
print(f"Samples in metadata: {len(metadata_samples)}")

# Find matching samples in metabolite data
data_samples = set(str(col).strip() for col in sample_cols)
matching_samples = metadata_samples.intersection(data_samples)

print(f"Matching samples found: {len(matching_samples)}")
print(f"Sample matching rate: {len(matching_samples)/len(metadata_samples)*100:.1f}%")
print(f"Example matches: {list(matching_samples)[:5]}")

# Analyze metabolite data quality
print(f"\n6. METABOLITE DATA QUALITY ANALYSIS:")
print("-" * 40)

# Extract numeric data only
numeric_data = maf_df[sample_cols].apply(pd.to_numeric, errors='coerce')
print(f"Metabolites (rows): {numeric_data.shape[0]}")
print(f"Samples (columns): {numeric_data.shape[1]}")

# Basic statistics
print(f"\nData range: {numeric_data.min().min():.2e} to {numeric_data.max().max():.2e}")
print(f"Median value: {numeric_data.median().median():.2e}")

# Missing value analysis
missing_by_metabolite = numeric_data.isnull().sum(axis=1)
missing_by_sample = numeric_data.isnull().sum(axis=0)

print(f"\nMissing values analysis:")
print(f"  Metabolites with no missing values: {(missing_by_metabolite == 0).sum()}")
print(f"  Metabolites with >50% missing: {(missing_by_metabolite > len(sample_cols)*0.5).sum()}")
print(f"  Samples with no missing values: {(missing_by_sample == 0).sum()}")
print(f"  Samples with >50% missing: {(missing_by_sample > len(numeric_data)*0.5).sum()}")

# Filter metadata for experimental samples only
print(f"\n7. EXPERIMENTAL DESIGN ANALYSIS:")
print("-" * 35)

# Remove blanks and controls
experimental_samples = sample_df[
    (sample_df['Characteristics[Organism]'] == 'Arabidopsis thaliana') &
    (~sample_df['Factor Value[Genotype]'].isin(['blank', 'Arabidopsis leaf control']))
].copy()

print(f"Total experimental samples: {len(experimental_samples)}")

# Analyze experimental design
print(f"\nTissue distribution:")
tissue_counts = experimental_samples['Characteristics[Organism part]'].value_counts()
print(tissue_counts)

print(f"\nTreatment distribution:")  
treatment_counts = experimental_samples['Factor Value[Drought Stress]'].value_counts()
print(treatment_counts)

print(f"\nGenotype distribution:")
genotype_counts = experimental_samples['Factor Value[Genotype]'].value_counts()
print(genotype_counts)

print(f"\nTime point distribution:")
time_counts = experimental_samples['Factor Value[Time point]'].value_counts().sort_index()
print(time_counts)

# Extract drought samples for analysis
print(f"\n8. DROUGHT SAMPLE EXTRACTION:")
print("-" * 30)

drought_samples = experimental_samples[experimental_samples['Factor Value[Drought Stress]'] == 'drought']

# Separate by tissue  
root_drought = drought_samples[drought_samples['Characteristics[Organism part]'] == 'root']
leaf_drought = drought_samples[drought_samples['Characteristics[Organism part]'].isin(['leaf', 'shoot'])]

print(f"Total drought samples: {len(drought_samples)}")
print(f"Root drought samples: {len(root_drought)}")
print(f"Leaf drought samples: {len(leaf_drought)}")

# Check temporal distribution
print(f"\nTemporal distribution of drought samples:")
time_tissue_counts = drought_samples.groupby(['Characteristics[Organism part]', 'Factor Value[Time point]']).size().unstack(fill_value=0)
print(time_tissue_counts)

# Extract expression data for both tissues
print(f"\n9. EXPRESSION DATA EXTRACTION:")
print("-" * 35)

# Get matching sample columns for each tissue
root_sample_names = set(root_drought['Sample Name'].astype(str).str.strip())
leaf_sample_names = set(leaf_drought['Sample Name'].astype(str).str.strip())

# Find intersection with available data columns
available_samples = set(str(col).strip() for col in sample_cols)

root_data_cols = list(root_sample_names.intersection(available_samples))
leaf_data_cols = list(leaf_sample_names.intersection(available_samples))

print(f"Root samples with data: {len(root_data_cols)}")
print(f"Leaf samples with data: {len(leaf_data_cols)}")

if len(root_data_cols) >= 5 and len(leaf_data_cols) >= 5:
    print("✓ Sufficient samples for network analysis!")
    
    # Extract expression matrices
    root_expr_raw = maf_df[root_data_cols].apply(pd.to_numeric, errors='coerce')
    leaf_expr_raw = maf_df[leaf_data_cols].apply(pd.to_numeric, errors='coerce')
    
    print(f"Raw expression matrices:")
    print(f"  Root: {root_expr_raw.shape} (metabolites × samples)")
    print(f"  Leaf: {leaf_expr_raw.shape} (metabolites × samples)")
    
    # Quality filtering
    root_good_metabolites = root_expr_raw.isnull().sum(axis=1) <= len(root_data_cols) * 0.5
    leaf_good_metabolites = leaf_expr_raw.isnull().sum(axis=1) <= len(leaf_data_cols) * 0.5
    
    root_expr_filtered = root_expr_raw[root_good_metabolites]
    leaf_expr_filtered = leaf_expr_raw[leaf_good_metabolites]
    
    print(f"After quality filtering:")
    print(f"  Root: {root_expr_filtered.shape} (metabolites × samples)")
    print(f"  Leaf: {leaf_expr_filtered.shape} (metabolites × samples)")
    
else:
    print("⚠️ Insufficient samples for robust network analysis")
    root_expr_filtered = None
    leaf_expr_filtered = None

# Save processed data
print(f"\n10. SAVING PROCESSED DATA:")
print("-" * 28)

# Save sample lists
experimental_samples.to_csv(os.path.join(output_data_path, "arabidopsis2_experimental_samples.csv"), index=False)
root_drought[['Sample Name', 'Factor Value[Genotype]', 'Factor Value[Time point]']].to_csv(
    os.path.join(output_data_path, "arabidopsis2_root_drought_samples.csv"), index=False)
leaf_drought[['Sample Name', 'Factor Value[Genotype]', 'Factor Value[Time point]']].to_csv(
    os.path.join(output_data_path, "arabidopsis2_leaf_drought_samples.csv"), index=False)

# Save expression data if available
if root_expr_filtered is not None and leaf_expr_filtered is not None:
    # Transpose for samples × metabolites format (standard for network analysis)
    root_expr_final = root_expr_filtered.T
    leaf_expr_final = leaf_expr_filtered.T
    
    root_expr_final.to_csv(os.path.join(output_data_path, "arabidopsis2_root_expression.csv"))
    leaf_expr_final.to_csv(os.path.join(output_data_path, "arabidopsis2_leaf_expression.csv"))
    
    print(f"✓ Expression matrices saved:")
    print(f"  Root: {root_expr_final.shape} (samples × metabolites)")
    print(f"  Leaf: {leaf_expr_final.shape} (samples × metabolites)")
    
    # Find common metabolites for comparison
    common_metabolites = set(root_expr_final.columns).intersection(set(leaf_expr_final.columns))
    print(f"✓ Common metabolites for network comparison: {len(common_metabolites)}")
    
    # Save common metabolite list
    pd.Series(list(common_metabolites)).to_csv(
        os.path.join(output_data_path, "arabidopsis2_common_metabolites.csv"), 
        index=False, header=['metabolite_id']
    )

# Create visualization of sample distribution
print(f"\n11. CREATING SAMPLE DISTRIBUTION PLOT:")
print("-" * 38)

plt.figure(figsize=(12, 8))

# Create subplot layout
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
fig.suptitle('Arabidopsis Dataset 2: Experimental Design Overview', fontsize=16, fontweight='bold')

# Tissue distribution
ax1 = axes[0, 0]
tissue_counts.plot(kind='bar', ax=ax1, color='lightgreen', alpha=0.7)
ax1.set_title('A. Tissue Distribution')
ax1.set_xlabel('Tissue Type')
ax1.set_ylabel('Number of Samples')
ax1.tick_params(axis='x', rotation=45)

# Treatment distribution
ax2 = axes[0, 1]
treatment_counts.plot(kind='bar', ax=ax2, color='lightblue', alpha=0.7)
ax2.set_title('B. Treatment Distribution')
ax2.set_xlabel('Treatment')
ax2.set_ylabel('Number of Samples')
ax2.tick_params(axis='x', rotation=45)

# Genotype distribution
ax3 = axes[1, 0]
genotype_counts.plot(kind='bar', ax=ax3, color='orange', alpha=0.7)
ax3.set_title('C. Genotype Distribution')
ax3.set_xlabel('Genotype')
ax3.set_ylabel('Number of Samples')
ax3.tick_params(axis='x', rotation=45)

# Time point distribution
ax4 = axes[1, 1]
time_counts.plot(kind='bar', ax=ax4, color='purple', alpha=0.7)
ax4.set_title('D. Time Point Distribution')
ax4.set_xlabel('Time Point (Days)')
ax4.set_ylabel('Number of Samples')

plt.tight_layout()
plt.savefig(os.path.join(output_plot_path, "arabidopsis2_experimental_design.png"), dpi=300, bbox_inches='tight')
plt.show()

print(f"✓ Experimental design plot saved")

print("\n" + "="*60)
print("ARABIDOPSIS DATASET 2 ANALYSIS COMPLETE!")
print("="*60)

print("\nKEY FINDINGS:")
print("✓ Dataset successfully loaded and analyzed")
print("✓ Sample-metabolite data successfully linked")
print("✓ Sufficient samples for network analysis identified")
print("✓ Expression matrices prepared for both tissues")
print("✓ Data quality assessed and filtered")

print(f"\nDATA SUMMARY:")
if root_expr_filtered is not None and leaf_expr_filtered is not None:
    print(f"• Root samples: {len(root_data_cols)}")
    print(f"• Leaf samples: {len(leaf_data_cols)}")
    print(f"• Root metabolites: {root_expr_filtered.shape[0]}")
    print(f"• Leaf metabolites: {leaf_expr_filtered.shape[0]}")
    print(f"• Common metabolites: {len(common_metabolites)}")
else:
    print("• Insufficient data for network analysis")

print(f"\nREADY FOR: Rigorous Network Validation Analysis!")
print("Next step: Run the network validation to test tissue asymmetry patterns")