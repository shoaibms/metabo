"""
ARABIDOPSIS DATASET 2 EXPLORER
==============================
Comprehensive exploration of the second Arabidopsis dataset to understand:
- File contents and formats
- Sample metadata (genotypes, tissues, timepoints)
- Metabolite data structure and quality
- Missing values and data distributions
- Experimental design suitability for validation

Adapted for new dataset location: arabidopsis2/
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set up paths for NEW dataset
base_path = r"C:\Users\USER\Desktop\data_chem_3_10"
data_path = os.path.join(base_path, "data", "v-data", "arabidopsis2")  # NEW LOCATION
output_plot_path = os.path.join(base_path, "validation", "output", "plot")
output_data_path = os.path.join(base_path, "validation", "output", "data")

# Create output directories if they don't exist
os.makedirs(output_plot_path, exist_ok=True)
os.makedirs(output_data_path, exist_ok=True)

print("="*60)
print("ARABIDOPSIS DATASET 2 EXPLORATION")
print("="*60)

# STEP 1: List all available files in NEW dataset
print("\n1. FILES IN NEW DATASET DIRECTORY:")
print("-" * 35)
try:
    files = os.listdir(data_path)
    for i, file in enumerate(files, 1):
        file_size = os.path.getsize(os.path.join(data_path, file)) / 1024  # KB
        print(f"{i:2d}. {file:<50} ({file_size:.1f} KB)")
    
    if not files:
        print("⚠️  No files found in the directory!")
        print(f"Please check if files exist in: {data_path}")
        exit()
        
except Exception as e:
    print(f"❌ Error accessing data directory: {e}")
    print(f"Please check if path exists: {data_path}")
    exit()

# STEP 2: Identify key files (flexible naming)
print("\n2. IDENTIFYING KEY METADATA FILES:")
print("-" * 35)

# Find sample metadata file (various possible names)
sample_files = []
assay_files = []
metabolite_files = []
investigation_files = []

for file in files:
    file_lower = file.lower()
    
    # Sample metadata files
    if file.startswith('s_') and file.endswith('.txt'):
        sample_files.append(file)
    elif 'sample' in file_lower and file.endswith(('.txt', '.csv', '.tsv')):
        sample_files.append(file)
    
    # Assay files
    if file.startswith('a_') and file.endswith('.txt'):
        assay_files.append(file)
    elif 'assay' in file_lower and file.endswith(('.txt', '.csv', '.tsv')):
        assay_files.append(file)
    
    # Metabolite data files
    if file.startswith('m_') and file.endswith(('.txt', '.tsv')):
        metabolite_files.append(file)
    elif any(keyword in file_lower for keyword in ['metabolite', 'maf', 'data', 'matrix']) and file.endswith(('.txt', '.csv', '.tsv')):
        metabolite_files.append(file)
    
    # Investigation files
    if file.startswith('i_') and file.endswith('.txt'):
        investigation_files.append(file)
    elif 'investigation' in file_lower and file.endswith(('.txt', '.csv', '.tsv')):
        investigation_files.append(file)

print(f"Sample metadata files: {sample_files}")
print(f"Assay files: {assay_files}")  
print(f"Metabolite files: {metabolite_files}")
print(f"Investigation files: {investigation_files}")

# STEP 3: Load and explore sample metadata (try all candidate files)
sample_df = None
for sample_file in sample_files:
    if sample_df is not None:
        break
        
    print(f"\n3. EXPLORING SAMPLE METADATA: {sample_file}")
    print("-" * 45)
    
    try:
        # Try different separators
        sample_df = pd.read_csv(os.path.join(data_path, sample_file), sep='\t', low_memory=False)
        print(f"✓ Successfully loaded with tab separator")
        break
    except:
        try:
            sample_df = pd.read_csv(os.path.join(data_path, sample_file), sep=',', low_memory=False)
            print(f"✓ Successfully loaded with comma separator")
            break
        except Exception as e:
            print(f"✗ Error loading sample file {sample_file}: {e}")
            continue

if sample_df is not None:
    print(f"Sample file shape: {sample_df.shape}")
    print(f"Columns: {list(sample_df.columns)}")
    print("\nFirst 3 rows:")
    print(sample_df.head(3))
    
    # Look for key columns (flexible matching)
    key_patterns = ['genotype', 'tissue', 'time', 'treatment', 'sample', 'condition', 'stress', 'drought']
    found_columns = []
    for col in sample_df.columns:
        for pattern in key_patterns:
            if pattern.lower() in col.lower():
                found_columns.append(col)
                break
    
    print(f"\nPotential key columns found: {found_columns}")
else:
    print("\n❌ Could not load any sample metadata file!")

# STEP 4: Explore assay metadata (if available)
if assay_files and sample_df is not None:
    assay_file = assay_files[0]
    print(f"\n4. EXPLORING ASSAY METADATA: {assay_file}")
    print("-" * 40)
    
    try:
        assay_df = pd.read_csv(os.path.join(data_path, assay_file), sep='\t', low_memory=False)
        print(f"✓ Successfully loaded assay file")
        print(f"Assay file shape: {assay_df.shape}")
        print(f"Columns: {list(assay_df.columns)}")
        print("\nFirst 2 rows:")
        print(assay_df.head(2))
    except Exception as e:
        print(f"✗ Error loading assay file: {e}")

# STEP 5: Look for metabolite data files (comprehensive search)
print(f"\n5. SEARCHING FOR METABOLITE DATA FILES:")
print("-" * 40)

# Prioritize files by likely importance
all_candidates = []
for file in files:
    if file.endswith(('.txt', '.csv', '.tsv', '.xlsx')):
        file_size = os.path.getsize(os.path.join(data_path, file))
        all_candidates.append((file, file_size))

# Sort by size (larger files more likely to contain data)
all_candidates.sort(key=lambda x: x[1], reverse=True)

print("All potential data files (sorted by size):")
for i, (file, size) in enumerate(all_candidates[:10], 1):  # Show top 10
    print(f"{i}. {file} ({size/1024:.1f} KB)")

# STEP 6: Try to load the most promising metabolite files
metabolite_df = None
loaded_file = None

# Try metabolite-specific files first, then largest files
candidate_files = metabolite_files + [f[0] for f in all_candidates[:5]]

for candidate_file in candidate_files:
    if metabolite_df is not None:
        break
        
    print(f"\n6. EXPLORING METABOLITE DATA: {candidate_file}")
    print("-" * 45)
    
    try:
        # Try loading with different parameters
        metabolite_df = pd.read_csv(os.path.join(data_path, candidate_file), sep='\t', low_memory=False)
        loaded_file = candidate_file
        print(f"✓ Successfully loaded metabolite file with tab separator")
        break
    except:
        try:
            metabolite_df = pd.read_csv(os.path.join(data_path, candidate_file), sep=',', low_memory=False)
            loaded_file = candidate_file
            print(f"✓ Successfully loaded metabolite file with comma separator")
            break
        except Exception as e:
            print(f"✗ Error loading {candidate_file}: {e}")
            continue

if metabolite_df is not None:
    print(f"Metabolite file shape: {metabolite_df.shape}")
    print(f"First 5 columns: {list(metabolite_df.columns[:5])}")
    print(f"Last 5 columns: {list(metabolite_df.columns[-5:])}")
    
    # Check if samples are in rows or columns
    print(f"\nFirst few values of first column:")
    print(metabolite_df.iloc[:5, 0])
    
    # Look for sample identifiers
    sample_like_cols = []
    for col in metabolite_df.columns[:10]:  # Check first 10 columns
        if any(keyword in str(col).lower() for keyword in ['sample', 'name', 'id']):
            sample_like_cols.append(col)
    print(f"Potential sample identifier columns: {sample_like_cols}")
    
    # Check data types and missing values
    numeric_cols = metabolite_df.select_dtypes(include=[np.number]).columns
    print(f"\nNumeric columns: {len(numeric_cols)}")
    
    if len(numeric_cols) > 0:
        missing_percent = (metabolite_df[numeric_cols].isnull().sum() / len(metabolite_df) * 100).describe()
        print("Missing values statistics (%):")
        print(missing_percent)

# STEP 7: DEEP DIVE DATA ANALYSIS
print(f"\n7. COMPREHENSIVE DATA ANALYSIS:")
print("-" * 40)

if sample_df is not None:
    print("\nSAMPLE METADATA DEEP DIVE:")
    print("-" * 30)
    
    # Detailed column analysis
    for col in sample_df.columns:
        unique_vals = sample_df[col].nunique()
        total_vals = len(sample_df[col])
        null_count = sample_df[col].isnull().sum()
        
        print(f"\nColumn: '{col}'")
        print(f"  Type: {sample_df[col].dtype}")
        print(f"  Unique values: {unique_vals}/{total_vals}")
        print(f"  Missing values: {null_count}")
        
        # Show unique values if categorical (< 20 unique values)
        if unique_vals < 20 and unique_vals > 0:
            unique_list = sample_df[col].dropna().unique()
            print(f"  Values: {list(unique_list)}")
        elif sample_df[col].dtype in ['object', 'string']:
            # Show first few examples for text columns
            examples = sample_df[col].dropna().head(3).tolist()
            print(f"  Examples: {examples}")

if metabolite_df is not None:
    print("\n\nMETABOLITE DATA DEEP DIVE:")
    print("-" * 30)
    
    # Data orientation check
    print(f"Data shape: {metabolite_df.shape} (rows × columns)")
    
    # Identify metadata vs data columns
    metadata_cols = []
    data_cols = []
    
    for col in metabolite_df.columns:
        if metabolite_df[col].dtype in ['object', 'string']:
            metadata_cols.append(col)
        elif pd.api.types.is_numeric_dtype(metabolite_df[col]):
            data_cols.append(col)
    
    print(f"Metadata columns: {len(metadata_cols)}")
    print(f"Numeric data columns: {len(data_cols)}")
    
    # Analyze metadata columns
    if metadata_cols:
        print(f"\nMETADATA COLUMNS ANALYSIS:")
        for col in metadata_cols[:5]:  # First 5 metadata columns
            unique_count = metabolite_df[col].nunique()
            print(f"  {col}: {unique_count} unique values")
            if unique_count < 10:
                print(f"    Values: {list(metabolite_df[col].dropna().unique())}")
            else:
                examples = metabolite_df[col].dropna().head(3).tolist()
                print(f"    Examples: {examples}")
    
    # Analyze numeric data
    if data_cols:
        print(f"\nNUMERIC DATA ANALYSIS:")
        numeric_data = metabolite_df[data_cols]
        
        # Overall statistics
        print(f"  Total metabolite features: {len(data_cols)}")
        print(f"  Data range: {numeric_data.min().min():.2e} to {numeric_data.max().max():.2e}")
        print(f"  Median value: {numeric_data.median().median():.2e}")
        
        # Missing value analysis
        missing_by_col = numeric_data.isnull().sum()
        print(f"\nMISSING VALUES ANALYSIS:")
        print(f"  Columns with no missing values: {(missing_by_col == 0).sum()}")
        print(f"  Columns with >50% missing: {(missing_by_col > len(numeric_data)*0.5).sum()}")
        print(f"  Average missing per column: {missing_by_col.mean():.1f}")
        
        # Data distribution analysis
        print(f"\nDATA DISTRIBUTION ANALYSIS:")
        non_zero_counts = (numeric_data > 0).sum()
        print(f"  Columns with >80% non-zero values: {(non_zero_counts > len(numeric_data)*0.8).sum()}")
        print(f"  Columns with <20% non-zero values: {(non_zero_counts < len(numeric_data)*0.2).sum()}")

# STEP 8: Sample identification strategy
print(f"\n8. SAMPLE IDENTIFICATION STRATEGY:")
print("-" * 40)

if sample_df is not None and metabolite_df is not None:
    print("Attempting to link sample metadata with metabolite data...")
    
    # Look for common identifiers
    sample_id_cols = [col for col in sample_df.columns if 'sample' in col.lower() or 'name' in col.lower()]
    metabolite_id_cols = [col for col in metabolite_df.columns if 'sample' in col.lower() or 'name' in col.lower()]
    
    print(f"Sample ID columns in metadata: {sample_id_cols}")
    print(f"Sample ID columns in metabolite data: {metabolite_id_cols}")
    
    # Check if sample names match between files
    if sample_id_cols and metabolite_id_cols:
        sample_names_meta = set(sample_df[sample_id_cols[0]].dropna().astype(str))
        sample_names_data = set(metabolite_df[metabolite_id_cols[0]].dropna().astype(str))
        
        overlap = sample_names_meta.intersection(sample_names_data)
        print(f"Overlapping sample identifiers: {len(overlap)}")
        if len(overlap) > 0:
            print(f"Example matches: {list(overlap)[:3]}")

# STEP 9: Tissue and treatment identification
print(f"\n9. TISSUE & TREATMENT IDENTIFICATION:")
print("-" * 40)

if sample_df is not None:
    # Look for tissue information
    tissue_cols = [col for col in sample_df.columns if any(keyword in col.lower() for keyword in ['tissue', 'organ', 'part'])]
    
    if tissue_cols:
        for col in tissue_cols:
            tissues = sample_df[col].dropna().unique()
            print(f"Tissue column '{col}': {list(tissues)}")
    
    # Look for treatment/condition information
    treatment_cols = [col for col in sample_df.columns if any(keyword in col.lower() for keyword in ['treatment', 'condition', 'stress', 'drought'])]
    
    if treatment_cols:
        for col in treatment_cols:
            treatments = sample_df[col].dropna().unique()
            print(f"Treatment column '{col}': {list(treatments)}")
    
    # Look for genotype information
    genotype_cols = [col for col in sample_df.columns if any(keyword in col.lower() for keyword in ['genotype', 'line', 'strain', 'variety'])]
    
    if genotype_cols:
        for col in genotype_cols:
            genotypes = sample_df[col].dropna().unique()
            print(f"Genotype column '{col}': {list(genotypes)}")
    
    # Look for time information
    time_cols = [col for col in sample_df.columns if any(keyword in col.lower() for keyword in ['time', 'day', 'hour', 'period'])]
    
    if time_cols:
        for col in time_cols:
            timepoints = sorted(sample_df[col].dropna().unique())
            print(f"Time column '{col}': {list(timepoints)}")

print(f"\n10. SUMMARY AND RECOMMENDATIONS:")
print("-" * 35)
print("✓ COMPREHENSIVE data exploration complete!")
print(f"✓ Output will be saved to: {output_plot_path}")
print(f"✓ Data files will be saved to: {output_data_path}")

print("\nKEY FINDINGS TO REVIEW:")
print("1. Sample metadata structure and key variables")
print("2. Metabolite data dimensions and orientation") 
print("3. Missing value patterns and data quality")
print("4. Tissue types, treatments, and genotypes available")
print("5. Sample identifier linking strategy")

print("\nNext steps:")
print("1. Review the detailed output above")
print("2. Confirm tissue separation (root vs shoot)")
print("3. Identify drought vs control conditions")  
print("4. Proceed with network analysis pipeline")

print("\n" + "="*60)
print("ARABIDOPSIS DATASET 2 EXPLORATION COMPLETE!")
print("="*60)